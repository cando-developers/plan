(in-package :metal-binder)

(defvar *force-field* (foldamer:load-force-field t))

(defun load-rotamers ()
  (let ((rotamers "~/work/spiros/data/rotamers.cando"))
    (format t "About to load ~s~%" rotamers)
    (time (defvar *rotamer-db* (cando.serialize:load-cando rotamers)))
    (format t "Loaded.~%")))

(defparameter *olig-space*
  (topology:make-oligomer-space
   spiros:*spiros*
   `((spiros:ring-aminoacid :label :ringc)
     :amide ((spiros:apro4ss spiros:apro4rr)) (:side spiros:bnz) (:sideamide spiros::sideamide-ace)
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side ((spiros:3pr spiros:4pr) :label :pr1))
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side spiros:bnz)
     :dkp spiros::ampross
     :amide spiros::aminoacid
     :amide ((spiros:apro4ss spiros:apro4rr)) (:side spiros:bnz) (:sideamide spiros:sideamide-ace)
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side ((spiros:3pr spiros:4pr) :label :pr2))
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side spiros:bnz)
     :dkp spiros::ampross
     :amide spiros::aminoacid
     :amide ((spiros:apro4ss spiros:apro4rr)) (:side spiros:bnz) (:sideamide spiros:sideamide-ace)
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side (spiros:bipy :label :bipy))
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side spiros:bnz)
     :dkp spiros::ampross ((ring :+amide :+ring :ringc))
     )
   ))


(defun add-restraints-to-energy-function (assembler)
  (let* ((oligomer (first (topology:oligomers assembler)))
         (oligomer-space (topology:oligomer-space oligomer))
         (mpy1 (first (gethash :pr1 (topology:labeled-monomers oligomer-space))))
         (mpy2 (first (gethash :pr2 (topology:labeled-monomers oligomer-space))))
         (mbipy (first (gethash :bipy (topology:labeled-monomers oligomer-space))))
         (pos-py1 (gethash mpy1 (topology:monomer-positions assembler)))
         (pos-py2 (gethash mpy2 (topology:monomer-positions assembler)))
         (pos-bipy (gethash mbipy (topology:monomer-positions assembler)))
         (aggregate (topology:aggregate assembler))
         (res-py1 (topology:at-position aggregate pos-py1))
         (res-py2 (topology:at-position aggregate pos-py2))
         (res-bipy (topology:at-position aggregate pos-bipy))
         (py1-n (chem:atom-with-name res-py1 :N))
         (py2-n (chem:atom-with-name res-py2 :N))
         (bipy-n3 (chem:atom-with-name res-bipy :N3))
         (bipy-n8 (chem:atom-with-name res-bipy :N8))
         (energy-function (topology:energy-function assembler))
         (stretch (chem:energy-function/get-stretch-component energy-function))
         (atomtable (chem:energy-function/atom-table energy-function))
         (force 1000.0)
         )
    (chem:energy-stretch/add-stretch-term stretch atomtable py1-n py2-n force 2.78)
    (chem:energy-stretch/add-stretch-term stretch atomtable py1-n bipy-n3 force 2.95)
    (chem:energy-stretch/add-stretch-term stretch atomtable py2-n bipy-n8 force 2.95)
    )
  )

(defclass mc-job (cando.serialize:serializable)
  ((team :initarg :team :accessor team)
   (node :initarg :node :accessor node)
   (oligomer-index :initarg :oligomer-index :accessor oligomer-index))
  )


(defun setup (&key (teams 1) (nodes-per-team 12) max-sequences verbose)
  (let* ((number-of-sequences (or max-sequences (topology:number-of-sequences *olig-space*)))
         (oligomer-index 0)
         (jobs nil))
    (loop named outer
          do (loop for team below teams
                   do (loop for node below nodes-per-team
                            for job = (make-instance 'mc-job
                                                     :oligomer-index oligomer-index
                                                     :team team
                                                     :node node)
                            do (push job jobs)
                            do (incf oligomer-index)
                            do (format t "oligomer-index team ~a node ~a ~a/~a~%" team node oligomer-index number-of-sequences)
                            when (= oligomer-index number-of-sequences)
                              do (return-from outer nil))))
    (loop for team below teams
          do (with-open-file (tfout (format nil "team~a.team" team) :direction :output)
               (loop for node below nodes-per-team
                     do (format tfout "$CLASP -l run.lisp -e \"(metal-binder:run ~a ~a)\"~%" team node))))
    (ensure-directories-exist "data/jobs.cando")
    (cando.serialize:save-cando jobs "data/jobs.cando")))


(defun run (team node &key verbose)
  (load-rotamers)
  (let* ((print-lock (bordeaux-threads:make-recursive-lock))
         (all-jobs (cando.serialize:load-cando "data/jobs.cando"))
         (jobs (loop for job in all-jobs
                     when (and (= (team job) team) (= (node job) node))
                       collect job))
         (file (pathname (format nil "~a/team-~anode-~a.output" "data" team node))))
    (with-open-file (fout file :direction :output)
      (format fout "(~%")
      (lparallel:pmapcar (lambda (job)
                           (let* ((oligomer-index job)
                                  (solution (do-monte-carlo *olig-space* oligomer-index :verbose verbose)))
                             (bordeaux-threads:with-recursive-lock-held (print-lock)
                               (let* ((*print-readably* t) 
                                      (*print-pretty* nil))
                                 (format fout "~s~%" solution))
                               (finish-output fout))))
                         jobs)
      (format fout ")~%")))
  )

(defclass mc-solution (cando.serialize:serializable)
  ((oligomer-index :initarg :oligomer-index :accessor oligomer-index)
   (score :initarg :score :accessor score)
   (rotamers :initarg :rotamers :accessor rotamers)))

(defmethod print-object ((obj mc-solution) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "oligomer-index ~a score: ~f" (oligomer-index obj) (score obj)))))

(defun do-monte-carlo (oligomer-space oligomer-index
                       &key (num-mc-runs 10)
                         (rotamer-db *rotamer-db*)
                         (verbose nil))
  (let* ((olig (topology:make-oligomer oligomer-space oligomer-index))
         (olig-shape (topology:make-oligomer-shape olig rotamer-db))
         (assembler (topology:make-assembler (list olig)))
         (energy-function (topology:energy-function assembler))
         (coords (topology:make-coordinates-for-assembler assembler))
         (bs (foldamer:make-backbone-rotamer-stepper olig-shape))
         (best-solution nil))
    (add-restraints-to-energy-function assembler)
    (loop for mc-index below num-mc-runs
          do (progn
               (topology:write-rotamers bs (foldamer:random-rotamers bs))
               (let* ((ss (foldamer:make-sidechain-rotamer-stepper olig-shape))
                      )
                 (topology:write-rotamers ss (foldamer:random-rotamers ss))
                 ;; Is this where I would add the distance restraints between amines?
                 (macrocycle:mopt-backbone olig-shape assembler coords)
                 (macrocycle:mopt-sidechain olig-shape assembler coords))
               (let* ((vec (topology:read-rotamers olig-shape))
                      (solution (make-instance 'mc-solution
                                               :oligomer-index oligomer-index
                                               :score (chem:evaluate-energy energy-function coords)
                                               :rotamers vec)))
                 (when (or (null best-solution) (< (score solution) (score best-solution)))
                   (when verbose (format t "Found a better solution ~s~%" solution))
                   (setf best-solution solution)))))
    best-solution
    )
  )


(defun analyze ()
  (let* ((files (directory "data/**.output"))
         (all-data (loop for file in files
                         append (with-open-file (fin file :direction :input)
                                  (let ((data (read fin)))
                                    data))))
         (sorted (sort all-data #'< :key #'score)))
    (with-open-file (fout "results.results" :direction :output)
      (format fout "(~%")
      (loop for solution in sorted
            do (let* ((*print-readably* t) 
                      (*print-pretty* nil))
                 (format fout "~s~%" solution)))
      (format fout ")~%")
      )))

(defun solution-aggregate (solution &key (rotamer-db *rotamer-db*))
  (let* ((olig (topology:make-oligomer *olig-space* (oligomer-index solution)))
         (olig-shape (topology:make-oligomer-shape olig rotamer-db))
         (assembler (topology:make-assembler (list olig)))
         (coords (topology:make-coordinates-for-assembler assembler))
         )
    (topology:write-rotamers olig-shape (rotamers solution))
    (topology:fill-internals-from-oligomer-shape-and-adjust assembler olig-shape)
    (topology:build-all-atom-tree-external-coordinates-and-adjust assembler coords)
    (topology::copy-joint-positions-into-atoms assembler coords)
    (topology:aggregate assembler)
    )
  )

(defparameter command nil)
(defparameter args nil)


(defun command-line-dispatch ()
  (let* ((raw-args  (loop for ii below (sys:argc) collect (sys:argv ii)))
         (script-pos (position "--script" raw-args :test #'string=)))
    (when script-pos
      (setf command (elt raw-args (+ 2 script-pos)))
      (setf args (mapcar #'read-from-string
                                 (nthcdr (+ script-pos 3)
                                         (loop for ii below (sys:argc) collect (sys:argv ii)))))))
  (cond
    ((string= command "args")
     (format t "args = ~s~%" (loop for ii below (sys:argc) collect (sys:argv ii)))
     (format t "args = ~s~%" args))
    ((string= command "setup")
     (apply 'setup args))
    ((string= command "run")
     (apply 'run args))
    ((string= command "analyze")
     (apply 'analyze args))
    (t (format t "Handle cmd ~s~%" command))
    )
  )


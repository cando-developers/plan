(in-package :metal-binder)

(defvar *num-mc-runs* 50)

(defvar *olig-space* nil)

(defun find-atom (assembler monomer-label atom-name)
  (let* ((oligomer-shape (first (topology:oligomer-shapes assembler)))
         (oligomer (topology:oligomer oligomer-shape))
         (oligomer-space (topology:oligomer-space oligomer))
         (monomer (gethash monomer-label (topology:labeled-monomers oligomer-space)))
         (pos (gethash monomer (topology:monomer-positions assembler)))
         (aggregate (topology:aggregate assembler))
         (residue (topology:at-position aggregate pos))
         (atom (chem:atom-with-name residue atom-name))
         )
    atom))

(defgeneric add-restraints-to-energy-function (assembler))

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


(defun run (team node &key verbose (parallel t))
  (load-rotamers)
  (let* ((print-lock (bordeaux-threads:make-recursive-lock))
         (all-jobs (cando.serialize:load-cando "data/jobs.cando"))
         (jobs (loop for job in all-jobs
                     when (and (= (team job) team) (= (node job) node))
                     collect job))
         (file (pathname (format nil "~a/team-~anode-~a.output" "data" team node))))
    (with-open-file (fout file :direction :output)
      (flet ((run-one (job)
               (let* ((oligomer-index (oligomer-index job))
                      (solution (do-monte-carlo *olig-space* oligomer-index :verbose verbose)))
                 (bordeaux-threads:with-recursive-lock-held (print-lock)
                   (let* ((*print-readably* t) 
                          (*print-pretty* nil))
                     (format fout "~s~%" solution))
                   (finish-output fout)))))
        (if parallel
            (lparallel:pmapcar #'run-one jobs)
            (mapcar #'run-one jobs))
        ))
    ))

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
                       &key (num-mc-runs *num-mc-runs*)
                         (rotamer-db *rotamer-db*)
                         (verbose t))
  (declare (special *rotamer-db*))
  (let* ((olig (topology:make-oligomer oligomer-space oligomer-index))
         (olig-shape (topology:make-oligomer-shape olig rotamer-db))
         (assembler (topology:make-assembler (list olig-shape)))
         (energy-function (topology:energy-function assembler))
         (coords (topology:make-coordinates-for-assembler assembler))
         (permissible-backbone-rotamers (topology:make-permissible-backbone-rotamers olig-shape))
         (best-solution nil))
    (add-restraints-to-energy-function assembler)
    (loop for mc-index below num-mc-runs
          do (format t "mc-index ~a oligomer-index: ~a~%" mc-index oligomer-index)
          do (loop
              (restart-case
                  (let ((rand-rots (topology:random-rotamers permissible-backbone-rotamers)))
                    (topology:write-rotamers olig-shape permissible-backbone-rotamers rand-rots)
                    (let* ((permissible-sidechain-rotamers (topology:make-permissible-sidechain-rotamers olig-shape)))
                      (topology:write-rotamers olig-shape permissible-sidechain-rotamers (topology:random-rotamers permissible-sidechain-rotamers))
                      ;; Restart the mopt 
                      (macrocycle:mopt-backbone olig-shape assembler coords :verbose verbose)
                      (macrocycle:mopt-sidechain olig-shape assembler coords :verbose verbose))
                    (return nil))
                (macrocycle:restart-monte-carlo ()
                  (format t "WARNING: metal-binder.lisp do-monte-carlo restart-monte-carlo was invoked~%"))))
          do (let* ((vec (topology:read-oligomer-shape-rotamers olig-shape))
                    (solution (make-instance 'mc-solution
                                             :oligomer-index oligomer-index
                                             :score (chem:evaluate-energy energy-function coords)
                                             :rotamers vec)))
               (when (or (null best-solution) (< (score solution) (score best-solution)))
                 (when verbose (format t "Found a better solution ~s~%" solution))
                 (setf best-solution solution))))
    best-solution
    )
  )


(defun analyze (&key (data "data") output)
  (let* ((files (directory (format nil "~a/*.output" data)))
         (all-data (let ((*readtable* cando.serialize::*cando-reader*))
                     (loop for file in files
                           do (progn
                                (format t "Reading ~s~%" file)
                                (finish-output t))
                           append (with-open-file (fin file :direction :input)
                                    (loop for one = (read fin nil :eof)
                                          when (eq one :eof)
                                            do (return results)
                                          collect one into results)))))
         (sorted (sort all-data #'< :key #'score)))
    (when output
      (with-open-file (fout output :direction :output)
        (loop for solution in sorted
              do (let* ((*print-readably* t) 
                        (*print-pretty* nil))
                   (format fout "~s~%" solution)))
        ))
    sorted
    ))

(defun solution-aggregate (solution)
  (let* ((rotamer-db (load-rotamers))
         (olig (topology:make-oligomer *olig-space* (oligomer-index solution)))
         (olig-shape (topology:make-oligomer-shape olig rotamer-db))
         (assembler (topology:make-assembler (list olig-shape)))
         (coords (topology:make-coordinates-for-assembler assembler))
         )
    (topology:write-oligomer-shape-rotamers olig-shape (rotamers solution))
    (topology:update-internals assembler olig-shape)
    (topology:build-all-atom-tree-external-coordinates-and-adjust assembler coords)
    (topology::copy-all-joint-positions-into-atoms assembler coords)
    (topology:aggregate assembler)))

(defun solution-sequence (solution)
  (let* ((olig (topology:make-oligomer *olig-space* (oligomer-index solution)))
         )
    (topology::canonical-sequence olig)))

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


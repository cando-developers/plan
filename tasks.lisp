(in-package :plan)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *test* nil)


(defclass plan-file (task:file)
  ((plan-name :initarg :plan-name :reader plan-name)))

(defclass search-file (plan-file)
  ((oligomer-index :initarg :oligomer-index :reader oligomer-index)
   (search-index :initarg :search-index :reader search-index)))

(defmethod task:file-pathname ((file search-file))
  (make-pathname :name (format nil "search-~D-~D" (oligomer-index file) (search-index file)) :type "cando"
                 :directory (list :relative "output")))

(defclass hit-file (plan-file)
  ((oligomer-index :initarg :oligomer-index :reader oligomer-index)
   (search-index :initarg :search-index :reader search-index)))

(defmethod task:file-pathname ((file hit-file))
  (make-pathname :name (format nil "hits-~D-~D" (oligomer-index file) (search-index file)) :type "cando"
                 :directory (list :relative "output")))

(defclass aggregate-file (plan-file)
  ())

(defmethod task:file-pathname ((file aggregate-file))
  (make-pathname :name "results" :type "cando"
                 :directory (list :relative "output")))

(defclass search-task (task:task)
  ((oligomer-space :initarg :oligomer-space :reader oligomer-space)
   (oligomer-index :initarg :oligomer-index :reader oligomer-index)
   (search-index :initarg :search-index :reader search-index)))

(defclass aggregate-task (task:task)
  ((oligomer-space :initarg :oligomer-space :reader oligomer-space)))

(defun build-task-graph (plan stream &key test)
  (format stream "build-task-graph in plan~%")
  (with-accessors ((name name)
                   (oligomer-space oligomer-space)
                   (search-count search-count))
      plan
    (let ((graph (make-instance 'task:task-graph)))
      (loop for oligomer-index below (if test
                                         1
                                         (topology:number-of-sequences oligomer-space))
            for oligomer = (topology:make-oligomer oligomer-space oligomer-index)
            do (loop for search-index below search-count
                     for input-file = (task:make-file graph 'search-file
                                                      :plan-name (name plan)
                                                      :oligomer-index oligomer-index
                                                      :search-index search-index)
                     for input-file-pathname = (task:file-pathname input-file)
                     for output-files = (task:make-files graph `(hit-file
                                                                 :plan-name ,(name plan)
                                                                 :oligomer-index ,oligomer-index
                                                                 :search-index ,search-index))
                     do (task:ensure-input-file-exists input-file)
                     do (task:make-task graph 'search-task
                                        :oligomer-space oligomer-space
                                        :oligomer-index oligomer-index
                                        :inputs (list input-file)
                                        :outputs output-files)))
      ;; Now a task to aggregate everything
      (let ((input-files (task:find-files graph :test (lambda (file)
                                                        (typep file 'hit-file))))
            (output-file (task:make-file graph 'aggregate-file
                                         :plan-name (name plan))))
        (task:make-task graph 'aggregate-task
                        :oligomer-space oligomer-space
                        :inputs input-files
                        :outputs (list output-file))
        graph))))

(defmethod task:execute ((task search-task))
  (let* ((output (first (task:outputs task)))
         (oligomer-space (oligomer-space task))
         (oligomer-index (oligomer-index task))
         (best-hit (if *test*
                       (do-monte-carlo oligomer-space oligomer-index :verbose t
                         :num-mc-runs 1)
                       (do-monte-carlo oligomer-space oligomer-index :verbose t))))
    ;; For now only do one search and skip multiple searches
    (loop for search below 0
          for hit = (do-monte-carlo oligomer-space oligomer-index :verbose t)
          when (< (score hit) (score best-hit))
            do (setf best-hit hit))
    (cando.serialize:save-cando best-hit (task:file-pathname output))))


(defmethod task:execute ((task aggregate-task))
  (let* ((inputs (task:inputs task))
         (output (first (task:outputs task)))
         (oligomer-space (oligomer-space task)))
    (let ((hits (loop for input in inputs
                      for hit = (cando.serialize:load-cando (task:file-pathname input))
                      collect hit)))
      (let* ((sorted-hits (sort hits #'< :key #'score))
             (results (make-instance 'mc-solution-set
                                     :oligomer-space oligomer-space
                                     :solutions sorted-hits)))
        (cando.serialize:save-cando results (task:file-pathname output))))))

(defun make-connection-pathname (plan-name)
  (let ((pn (plan-pathname plan-name)))
    (make-pathname :name "connection" :type "config" :defaults pn)))

(defun start-server (&rest args &key (show-remaining-task-limit 100) (log-to-file t) threaded to-stage endpoint)
  (let* ((plan (load-plan))
         (task-graph-callback (lambda (stream)
                                (build-task-graph plan stream)))
         (connection-path (make-connection-pathname nil)))
    (apply 'task:make-server task-graph-callback :connection-path connection-path args)))

(defun define-scoring-method (plan)
  (let* ((scorer (first (scorers plan)))
         (args (second scorer))
         (body (cddr scorer)))
    (eval `(defmethod plan:add-restraints-to-energy-function ,args ,@body))))


(defun start-client (&rest args &key (log-to-file t) to-stage threaded dont-execute)
  (let* ((plan (load-plan))
         (task-graph-callback (lambda (stream)
                                (build-task-graph plan stream)))
         (connection-path (make-connection-pathname nil)))
    (define-scoring-method plan)
    (apply 'task:make-client task-graph-callback :connection-path connection-path args)))

(defun run (&rest args &key test)
  (let* ((*test* test)
         (plan (load-plan))
         (task-graph-callback (lambda (stream) (build-task-graph plan stream :test test))))
    (define-scoring-method plan)
    (apply 'task:run-tasks task-graph-callback args)))

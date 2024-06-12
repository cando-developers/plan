(in-package :metal-binder)

(defclass search-file (task:file)
  ((oligomer-index :initarg :oligomer-index :reader oligomer-index)
   (search-index :initarg :search-index :reader search-index)))

(defmethod task:file-pathname ((file search-file))
  (make-pathname :name (format nil "search-~D-~D" (oligomer-index file) (search-index file)) :type "cando"
                 :directory '(:relative "output")))


(defclass hit-file (task:file)
  ((oligomer-index :initarg :oligomer-index :reader oligomer-index)
   (search-index :initarg :search-index :reader search-index)))

(defmethod task:file-pathname ((file hit-file))
  (make-pathname :name (format nil "hits-~D-~D" (oligomer-index file) (search-index file)) :type "cando"
                 :directory '(:relative "output")))


(defclass aggregate-file (task:file)
  ())

(defmethod task:file-pathname ((file aggregate-file))
  #P"output/results.cando")

(defclass search-task (task:task)
  ((oligomer-space :initarg :oligomer-space :reader oligomer-space)
   (oligomer-index :initarg :oligomer-index :reader oligomer-index)
   (search-index :initarg :search-index :reader search-index)))

(defclass aggregate-task (task:task)
  ((oligomer-space :initarg :oligomer-space :reader oligomer-space)))

(defun build-task-graph (oligomer-space &optional (searches 1))
  (let ((graph (make-instance 'task:task-graph)))
    (loop for oligomer-index below (topology:number-of-sequences oligomer-space)
          for oligomer = (topology:make-oligomer oligomer-space oligomer-index)
          do (loop for search-index below searches
                   for input-file = (task:make-file graph `(search-file
                                                              :oligomer-index ,oligomer-index
                                                              :search-index ,search-index))
                   for input-file-pathname = (task:file-pathname input-file)
                   for output-files = (task:make-files graph `(hit-file
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
          (output-file (task:make-file graph `(aggregate-file))))
      (task:make-task graph 'aggregate-task
                      :inputs input-files
                      :outputs (list output-file))
    graph)))

(defmethod task:execute ((task search-task))
  (let* ((inputs (task:inputs task))
         (output (first (task:outputs task)))
         (search-file (first inputs))
         (oligomer-space (oligomer-space task))
         (oligomer-index (oligomer-index task))
         (best-hit (do-monte-carlo oligomer-space oligomer-index :verbose nil)))
    (loop for search below 0
          for hit = (do-monte-carlo oligomer-space oligomer-index :verbose nil)
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


(defun make-server (&rest args &key (show-remaining-task-limit 100) (log-to-file t) connection-path threaded to-stage endpoint)
  (let ((task-graph (build-task-graph metal-binder:*olig-space*)))
    (apply 'task:make-server task-graph args)))

(defun make-client (&rest args &key (log-to-file t) connection-path to-stage threaded dont-execute)
  (let ((task-graph (build-task-graph metal-binder:*olig-space*)))
    (apply 'task:make-client task-graph args)))

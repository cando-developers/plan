(in-package :asdf-user)

(defsystem "metal-binder"
  :description "Build metal-binding molecules"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:spiros :macrocycle :foldamer)
  :serial t
  :components ((:file "packages")
               (:file "metal-binder")
               (:file "tasks")
               ))

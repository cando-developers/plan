(in-package :asdf-user)

(defsystem "plan-demo"
  :description "Build metal-binding molecules"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:plan)
  :serial t
  :components ((:file "packages-demo")
               (:file "plan-demo")
               ))

(in-package :asdf-user)

(defsystem "metal-binder-demo"
  :description "Build metal-binding molecules"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:metal-binder)
  :serial t
  :components ((:file "packages-demo")
               (:file "metal-binder-demo")
               ))

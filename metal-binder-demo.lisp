(in-package :metal-binder-demo)

(setq metal-binder:*olig-space*
  (topology:make-oligomer-space
   spiros:*spiros*
   `((spiros:ring-aminoacid :label :ringc)
     :amide ((spiros:apro4ss spiros:apro4rr)) (:side spiros:bnz) (:amide spiros::ace)
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side ((spiros:3pr spiros:4pr) :label :pr1))
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side spiros:bnz)
     :dkp spiros::ampross
     :amide spiros::aminoacid
     :amide ((spiros:apro4ss spiros:apro4rr)) (:side spiros:bnz) (:amide spiros:ace)
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side ((spiros:3pr spiros:4pr) :label :pr2))
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side spiros:bnz)
     :dkp spiros::ampross
     :amide spiros::aminoacid
     :amide ((spiros:apro4ss spiros:apro4rr)) (:side spiros:bnz) (:amide spiros:ace)
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side (spiros:bipy :label :bipy))
     :dkp ((spiros:pro4ss spiros:pro4rr)) (:side spiros:bnz)
     :dkp spiros::ampross ((ring :+amide :+rev-amide :ringc))
     )
   ))


`((spiros:))

(defmethod metal-binder:add-restraints-to-energy-function (assembler)
  (let* ((py1-n (metal-binder:find-atom assembler :pr1 :N))
         (py2-n (metal-binder:find-atom assembler :pr2 :N))
         (bipy-n3 (metal-binder:find-atom assembler :bipy :N3))
         (bipy-n8 (metal-binder:find-atom assembler :bipy :N8))
         (energy-function (topology:energy-function assembler))
         (atomtable (chem:energy-function/atom-table energy-function))
         (stretch (chem:energy-function/get-stretch-component energy-function))
         (force 1000.0)
         )
    (chem:energy-stretch/add-stretch-term stretch atomtable py1-n py2-n force 2.78)
    (chem:energy-stretch/add-stretch-term stretch atomtable py1-n bipy-n3 force 2.95)
    (chem:energy-stretch/add-stretch-term stretch atomtable py2-n bipy-n8 force 2.95)
    )
  )

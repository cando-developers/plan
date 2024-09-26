(in-package :plan-demo)

(progn
  (defparameter bbo '(spiros:pro4ss spiros:pro4rr))
  (defparameter bbe '(spiros:pre4ss spiros:pre4rr)))

#+(or)
(progn
  (defparameter bbo '(spiros:pro4ss))
  (defparameter bbe '(spiros:pre4ss)))

(setq plan:*olig-space*
      (topology:make-oligomer-space
       spiros:*spiros*
       `((spiros:root-adabs :label :ringc)
         :dkp (,bbo) (:side spiros:bnz)
         :dkp (,bbo) (:side ((spiros:3pr spiros:4pr) :label :pr1))
         :dkp (,bbe) (:eside spiros:ebnz)
         :rev-amide spiros:rev-bala
         :rev-amide spiros:adabs
         :dkp (,bbo) (:side spiros:bnz)
         :dkp (,bbo) (:side ((spiros:3pr spiros:4pr) :label :pr2))
         :dkp (,bbe) (:eside spiros:ebnz)
         :rev-amide spiros:rev-bala
         :rev-amide spiros:adabs
         :dkp (,bbo) (:side spiros:bnz)
         :dkp (,bbo) (:side ((spiros:bipy) :label :bipy))
         :dkp (,bbe) (:eside spiros:ebnz)
         :rev-amide spiros:rev-bala ((ring :+rev-amide :+amide :ringc))
       )))

(defmethod plan:add-restraints-to-energy-function (assembler)
  (let* ((py1-n (plan:find-atom assembler :pr1 :N))
         (py2-n (plan:find-atom assembler :pr2 :N))
         (bipy-n3 (plan:find-atom assembler :bipy :N3))
         (bipy-n8 (plan:find-atom assembler :bipy :N8))
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

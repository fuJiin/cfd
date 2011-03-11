(ns cfd.odppde
  (:use [cfd.core]))

;; Explicit Methods  
(defn ftcs
  "Uses forward time/central space method to calculate n+1.
   Takes a 1x3 vector as grid."
  [grid x-step t-step alpha]
  (let [i-1     (first grid)
        i       (second grid)
        i+1     (last grid)
        coeff   (calc-coeff x-step t-step alpha)]
    ($= i + coeff * (i+1 - (2 * i) + i-1))))

(defn calc-coeff
  "Calculates common coefficent term used in explicit methods"
  [x-step t-step alpha])
    
(defn ftcs-stable?
  "Check ftcs stability condition"
  [x-step t-step alpha]
  (<= (calc-coeff x-step t-step alpha)
      1/2)))
      
(defn dufort-frankel
  "Uses DuFort-Frankel method to calculate n+1
   Takes a 2x3 matrix (2 1x3 vectors) as grid."
   [grid x-step t-step alpha]
   (let [n-1    (second (last grid))
         i-1    (first (first grid))
         i+1    (last (first grid))
         coeff  (calc-coeff x-step t-step alpha)]
    ($= ((1 - 2 * coeff) * n-1 + (2 * coeff) * (i+1 + i-1))
        /
        (1 + 2 * coeff))))

;; Implicit Methods
(defn laasonen
  "Uses Laasonen implicit method to solve n+1"
  [grid x-step t-step alpha])
  
(defn crank-nicolson
  "Uses Crank-Nicolson to solve n+1"
  [grid x-step t-step alpha])
  
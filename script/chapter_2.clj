(ns cfd.script.chapter-2
  (:use [cfd.fd-approx])
  (:use [incanter.core]))

(defn f [x]
  ($= x ** 3 - 5 * x))
  
(defn f-deriv [x]
  ($= 3 * (x ** 2) - 5))

(defn f-deriv-2 [x]
  ($= 6 * x))
  
(def meths [:backward :central :forward])
  
(def i-vals [0.5 1.5])  
  
(def step-sizes [0.00001 0.0001 0.001 0.01 0.1 0.2 0.3])

(defn to-table
  "Creates a table for specified finite difference approximation"
  [power order func i steps]
  (let [c-map   (get-map order)
        data    (map (fn [method]
                       (map #(fd-approx power order method func i %) steps))
                     (vec (keys c-map)))])
  
  (let [approx  (map #(fd-approx power order method func i %) steps)
        deriv   (case   power
                 1      f-deriv
                 2      f-deriv-2
                 :else  (throw (Exception. "Invalid power")))
        actual  (deriv i)
        headers (conj step-sizes "actual")
        data    (conj (vec approx) actual)]
    (dataset headers [data])))
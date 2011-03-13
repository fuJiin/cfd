(ns cfd.script.chapter-2
  (:use [cfd.fd-approx])
  (:use [incanter.core]))

(defn f [x]
  ($= x ** 3 - 5 * x))
  
(defn f-deriv [x]
  ($= 3 * (x ** 2) - 5))

(defn f-deriv-2 [x]
  ($= 6 * x))
  
(def i-vals [0.5 1.5])  
  
(def step-sizes [0.00001 0.0001 0.001 0.01 0.1 0.2 0.3])

(defn get-approx
  "Creates a table for specified finite difference approximation.
   Takes in power, order, function, x, and steps to iterate through.
   It'll find all methods available for the order, and format results
   for each method in a separate row."
  [power order func x steps]
  (let [c-map   (get-coeff-map order)
        data    (map (fn [method]
                       (cons (name method)
                             (map #(fd-approx power order method func x %) steps)))
                     (keys c-map))
        header  (cons "Differencing" step-sizes)]
    (dataset header data)))
    
(defn approx-compare
  "Runs approximations and compare to actual values"
  [power order func x steps]
  (let [approx  (get-approx power order func x steps)
        deriv   (case power
                  1   f-deriv
                  2   f-deriv-2)
        actual  (if deriv (deriv x) nil)]
    {:actual actual
     :approx approx}))
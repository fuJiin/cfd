(ns cfd.fd-approx
  (:use [incanter.core]))
  
;; Coefficients  
(def first-order-coeffs
  {:forward   [{0   -1
                1   1}
                
               {0   1
                1   -2
                2   1}
                
               {0   -1
                1   3
                2   -3
                3   1}
                
               {0   1
                1   -4
                2   6
                3   -4
                4   1}]
                
   :backward  [{-1  -1
                0   1}
                
               {-2  -1
                -1  -2
                0   1}
                
               {-3  -1
                -2  3
                -1  -3
                0   1}
                
               {-4  1
                -3  -4
                -2  6
                -1  -4
                0   1}]})
                
(def second-order-coeffs
  {:central   [{:denom  2
                -1      -1
                0       0
                1       1}
                
               {-1      1
                0       -2
                1       1}
                
               {:denom  2
                -2      -1
                -1      2
                0       0
                1       -2
                2       1}
                
               {-2      1
                -1      -4
                0       6
                1       -4
                2       1}]
    
   :forward   [{:denom  2
                0       -3
                1       4
                2       -1}
                
               {0       2
                1       -5
                2       4
                3       -1}
                
               {:denom  2
                0       -5
                1       18
                2       -24
                3       14
                4       -3}
                
               {0       3
                1       -14
                2       26
                3       -24
                4       11
                5       -2}]
                
   :backward  [{:denom  2
                -2      1
                -1      -4
                0       3}
                
               {-3      -1
                -2      4
                -1      -5
                0       2}
                
               {:denom  2
                -4      3
                -3      -14
                -2      24
                -1      -18
                0       5}
                
               {-5      -2
                -4      11
                -3      -24
                -2      26
                -1      -14
                0       3}]})
                
(def fourth-order-coeffs
  {:central   [{:denom  12
                -2      1
                -1      -8
                0       0
                1       8
                2       -1}
                
               {:denom  12
                -2      -1
                -1      16
                0       -30
                1       16
                2       -1}
                
               {:denom  8
                -3      1
                -2      -8
                -1      13
                0       0
                1       -13
                2       8
                3       -1}
                
               {:denom  6
                -3      -1
                -2      12
                -1      -39
                0       56
                1       -39
                2       12
                3       -1}]})

;; Approximation functions
(defn get-coeff-map 
  "Pick coefficients based on order"
  [order]
  (case order
    1   first-order-coeffs
    2   second-order-coeffs
    4   fourth-order-coeffs
    :else (throw (Exception. "Invalid approximation orde"))))

(defn derive-func
  "Returns a function that approximates the derivative of the
   specified order, for a given function.
   
   The returned function will take as arguments
   - func:  the function to derive derivative for
   - i:     position to calculate at
   - step:  grid steps"
  [power order method]
  (let [c-map   (get-map order)
        c-hash  (nth (c-map method) (dec power))  ;; narrow to specific method and power
        denom   (c-hash :denom)                   ;; filter out denominator
        coeffs  (dissoc c-hash :denom)]           ;; get actual coefficients
    (fn [func i step]
      (/
        (reduce +                                 ;; sum together...
          (map #(* (coeffs %)                     ;; a collection that multiplies coefficient...
                   (func (+ i (* % step))))       ;; ...with the value of function called on i+?
                (keys coeffs)))                   
        ($= step ** power)))))                    ;; divide by step expt power
                                                  
(defn fd-approx
  "Calculate approximated derivative of given order, using specified method,
   for a given function, position, and step size"
  [power order method func i step]
  (let [fd-deriv (derive-func power order method)]
    (fd-deriv func i step)))
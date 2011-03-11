(ns cfd.fd-approx
  (:use [incanter.core]))
  
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
                0   1}]
                
   :central   [{:denom  2
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
                2       1}]}})
                
(def sec-order-coeffs
  {:forward   [{:denom  2
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
  
(defn deriv-func
  [])
  
(defn fd-approx)
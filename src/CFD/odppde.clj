(ns cfd.odppde
  (:use [incanter.core]))

;; Explicit Methods  
(defn calc-coeff
  "Calculates common coefficent term used in explicit methods"
  [x-step t-step alpha])
  
(defn ftcs
  "Uses forward time/central space method to calculate n+1.
   Takes a 1x3 vector as grid."
  [grid x-step t-step alpha]
  (let [i-1     (first grid)
        i       (second grid)
        i+1     (last grid)
        coeff   (calc-coeff x-step t-step alpha)]
    ($= i + coeff * (i+1 - (2 * i) + i-1))))

(defn ftcs-stable?
  "Check ftcs stability condition"
  [x-step t-step alpha]
  (<= (calc-coeff x-step t-step alpha)
      1/2))
      
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
;; For simplication purposes, boundary conditions are assumed constant at edges
(defn diag-offset
  "Creates a diagonal matrix with offset.
   Positive offset shifts right; negative offset shifts left"
  ([diag-grid offset]
   (let [length (count diag-grid)]
     (matrix
       (for [i (range 0 length)]    ;; for each row
         (for [x (range 0 length)]  ;; iterate through columns
           (if (= i (- x offset))   ;; if at index with index
               (nth diag-grid i)    ;; assign grid value
               0))))))              ;; otherwise assign 0
  ([diag-grid]
    (diag-offset diag-grid 0)))     ;; default offset to 0 if no provided

(defn coeff-matrix
  "Create coefficients matrix from a, b, c values.
   Assumes coefficients to be constant for at all i values."
  [length a b c]
  (let [a-grid (diag (repeat length a))
        b-grid (diag-offset (repeat length b) 1)
        c-grid (diag-offset (repeat length c) -1)]
    (plus a-grid b-grid c-grid)))
              
(defn solve-with-tridiag
  "Solves implicit formula using tridiagonal matrix formulation.
   Takes in coefficents a, b, c, and a 1xn vector as RHS values for D"
   [a b c d-grid]
   (let [length (count d-grid)
         coeffs (coeff-matrix length a b c)]))
         
(defn laasonen
  "Uses Laasonen implicit method to solve n+1.
   Takes in 1xn grid of points"
  [grid x-step t-step alpha]
  (let [coeff   ($= alpha * t-step / (x-step ** 2))
        b       (- ($= 2 * coeff + 1))
        d-grid  (map #(- %) grid)] ;; RHS is negative of current point for Laasonen
    (solve-with-tridiag coeff b coeff d-grid)))
  
(defn crank-nicolson
  "Uses Crank-Nicolson to solve n+1"
  [grid x-step t-step alpha]
  (let [coeff   ($= alpha / (2 * (x-step ** 2)))
        b       (- ($= 2 * coeff + ($= t-step ** -1)))
        d-grid  (map (fn [pt]                   ;; mapping RHS of equation
                        (let [i (.indexOf grid pt)]
                          (if (or (= pt 0) (= ($= pt (.length grid) - 1)))
                              pt                ;; just return the point if at boundary conditions
                              ($= (pt / t-step) ;; otherwise calculate based on i-1, i, and i+1
                                  +
                                  (coeff * ($= (nth grid (inc i)) - 2 * pt + (nth grid (dec i))) 
                                               / 
                                               ($= x-step ** 2))))))
                      grid)]                    ;; apply map function to current grid
    (solve-with-tridiag coeff b coeff d-grid)))
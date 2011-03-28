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
  (let [a-grid (diag-offset (repeat length a) -1)
        b-grid (diag (repeat length b))
        c-grid (diag-offset (repeat length c) 1)]
    (plus a-grid b-grid c-grid)))

(defn solve-with-tridiag
  "Solves implicit formula using tridiagonal matrix formulation.
   Takes in coefficents a, b, c, and a 1xn vector as RHS values for D"
   [a b c d-grid orig-grid]
   (let [d-size (count d-grid)
         rhs    (remove nil?                            ;; filter out nils at beginning or end
                  (for [i (range 0 d-size)]
                    (cond
                      (or (= i 0) (= i (dec d-size)))   ;; return nil if at beginning or end
                        nil
                      (= i 1)                           ;; return D_n - a * u_1 for i = 2
                        ($= (nth d-grid i) - a * (nth orig-grid 0))
                      (= i (- d-size 2))                ;; return D_n - c * u_-1 for i = IM1
                        ($= (nth d-grid i) - c * (nth orig-grid (dec d-size)))
                      :else                             ;; otherwise return point on grid
                        (nth d-grid i))))
         length (count rhs)
         coeffs (coeff-matrix length a b c)]            ;; get coefficient matrix
      (vec
        (flatten
          (vector (nth orig-grid 0)                     ;; add initial boundary value with...
                  (mmult (solve coeffs) rhs)            ;; multiply inverse of coeff-matrix with rhs to solve grid
                  (nth orig-grid (dec d-size)))))))     ;; and end boundary value

(defn laasonen
  "Uses Laasonen implicit method to solve n+1.
   Takes in 1xn grid of points"
  [grid x-step t-step alpha]
  (let [coeff   ($= alpha * t-step / (x-step ** 2))
        b       (- ($= 2 * coeff + 1))
        d-grid  (map #(- %) grid)] ;; RHS is negative of current point for Laasonen
    (solve-with-tridiag coeff b coeff d-grid grid)))

(defn crank-nicolson
  "Uses Crank-Nicolson to solve n+1"
  [grid x-step t-step alpha]
  (let [coeff   ($= alpha / (2 * (x-step ** 2)))
        b       (- ($= 2 * coeff + (1 / t-step)))
        d-grid  (for [i (range 0 (count grid))]
                  (let [pt (nth grid i)]
                    (if (or (= i 0) (= i (dec (count grid))))
                      pt
                      (- ($= (pt / t-step)
                             +
                             (coeff * (nth grid (inc i)) - 2 * pt + (nth grid (dec i)))
                                      /
                                      (x-step ** 2))))))]
    (solve-with-tridiag coeff b coeff d-grid grid)))
;; One Dimensional Parabolic Partial Differential Equations.
;; All following function assume Dirichlet Boundary Conditions.
(ns cfd.odppde
  (:use [incanter.core]))

;; Explicit Methods
;; For simplication purposes, boundary conditions are assumed constant at edges.
;; This should probably be a macro later, so it can be expanded at compile time,
;; and not recalculated every iteration
(defn calc-coeff
  "Calculates common coefficent term used in explicit methods"
  [x-step t-step alpha]
  ($= alpha * t-step / (x-step ** 2)))

(defn ftcs
  "Uses forward time/central space method to calculate n+1.
   Takes a 1xn vector as grid."
  [grid x-step t-step alpha]
  (let [coeff   (calc-coeff x-step t-step alpha)
        length  (count grid)]
    (for [i (range 0 length)]
      (if (or (= 0 i) (= (dec length) i))
          (nth grid i)                  ;; return initial val at boundary conditions
          (let [pt      (nth grid i)    ;; otherwise calculate using FTCS diff approx
                next-pt (nth grid (inc i))
                prev-pt (nth grid (dec i))]
            ($= pt + coeff * (next-pt - (2 * pt) + prev-pt)))))))

(defn ftcs-stable?
  "Check ftcs stability condition"
  [x-step t-step alpha]
  (<= (calc-coeff x-step t-step alpha)
      1/2))

(defn dufort-frankel
  "Uses DuFort-Frankel method to calculate n+1,
   taking a 2xn matrix as grid.
   The length of each row in the 2xn matrix must match up."
   [grid x-step t-step alpha]
   (let [coeff     (calc-coeff x-step t-step alpha)
         pres-grid (last grid)
         prev-grid (first grid)
         length    (count pres-grid)]
    (for [i (range 0 length)]                   ;; iterate through the grid
      (if (or (= i 0) (= i (dec length)))
        (nth pres-grid i)                       ;; return initial val at boundary conditions
        (let [pt      (nth pres-grid i)         ;; otherwise calculate using dufort-frankel diff approx
              next-pt (nth pres-grid (inc i))
              prev-pt (nth pres-grid (dec i))
              past-pt (nth prev-grid i)]
          ($= ((1 - 2 * coeff) * past-pt + (2 * coeff) * (next-pt + prev-pt))
              /
              (1 + 2 * coeff)))))))

;; Implicit Methods
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
         rhs    (remove nil?                          ;; filter out nils at beginning or end
                  (for [i (range 0 d-size)]
                    (cond
                      (or (= i 0) (= i (dec d-size))) ;; return nil if at beginning or end
                        nil
                      (= i 1)                         ;; return D_n - a * u_1 for i = 2
                        ($= (nth d-grid i) - a * (nth orig-grid 0))
                      (= i (- d-size 2))              ;; return D_n - c * u_-1 for i = IM1
                        ($= (nth d-grid i) - c * (nth orig-grid (dec d-size)))
                      :else                           ;; otherwise return point on grid
                        (nth d-grid i))))
         length (count rhs)
         coeffs (coeff-matrix length a b c)]          ;; get coefficient matrix
      (vec
        (flatten
          (vector (nth orig-grid 0)                   ;; add initial boundary value with...
                  (mmult (solve coeffs) rhs)          ;; multiply inverse of coeff-matrix with rhs to solve grid
                  (nth orig-grid (dec d-size)))))))   ;; and end boundary value

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
  (let [b       (- ($= 2 + 2 / (alpha * t-step)))
        d-grid  (for [i (range 0 (count grid))]       ;; iterate through grid indices to generate RHS
                  (let [pt  (nth grid i)]
                    (if (or (= i 0) (= i (dec (count grid))))
                      pt                              ;; just return value at boundary conds
                      ($= ((0 - (nth grid (inc i)))   ;; otherwise calculate based on n-values
                            + ((2 - (2 * (x-step ** 2) / (alpha * t-step))) * pt)
                            - (nth grid (dec i)))
                          /
                          (x-step ** 2)))))]
    (solve-with-tridiag 1 b 1 d-grid grid)))          ;; solve with tridiagonal sys of eqns

; Iteration helpers
; Helps iterate through x and t steps, accumulating previous grid values
(defn grid-run
  "Iterate over grid using given function.  Accumulates results.
   ++lookback++ indicates how many levels of n-results to use when calculating n+1.
   This value should be 1 except when used with DuFort-Frankel method."
  [func t-limit init-grid x-step t-step alpha lookback]
  (loop [output     (if (= lookback 1)              ;; initiate output based on lookback
                        (conj [] init-grid)         ;; creates vec of a 1xn vec if lookback = 1
                        init-grid)                  ;; otherwise use initial grid provided
         grid       init-grid
         t-counter  0]
    (let [result    (func grid x-step t-step alpha) ;; call function on grid
          new-grid  (conj output (vec result))      ;; log result
          t         (inc t-counter)]                ;; add t
      (if (>= t t-limit)
          new-grid                                  ;; return new-grid if at t-limit
          (recur new-grid                           ;; otherwise recalculate next step
                 ;; get grid for next iteration based on lookback
                 ;; should use single 1xn vec if lookback is 1,
                 ;; otherwise use collection of m 1xn vectors to represent mxn matrix
                 (if (= lookback 1)
                     (last new-grid)
                     (subvec new-grid (- (count new-grid) lookback)))
                 t)))))

(defn regular-run
  "Shortcut method to do a grid-run with single lookback"
  [func t-limit init-grid x-step t-step alpha]
  (grid-run func t-limit init-grid x-step t-step alpha 1))

(defn dufort-frankel-run
  "Grid iteration for DuFort-Frankel method.
   Takes in 1xn grid, and generates first n+1 using FTCS if stable,
   otherwise uses Crank-Nicolson method."
  [t-limit init-grid x-step t-step alpha]
  (let [stable    (ftcs-stable? x-step t-step alpha)
        fut-grid  (if stable
                      (ftcs           init-grid x-step t-step alpha)
                      (crank-nicolson init-grid x-step t-step alpha))
        grid      (conj [] init-grid (vec fut-grid))]
    (grid-run dufort-frankel t-limit grid x-step t-step alpha 2)))

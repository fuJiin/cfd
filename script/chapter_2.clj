(ns cfd.script.chapter-2
  (:use [cfd.fd-approx])
  (:use [incanter core io]))

;; Problem Conditions ;;
(defn f [x]
  ($= x ** 3 - 5 * x))

(defn f-deriv [x]   ; first derivative
  ($= 3 * (x ** 2) - 5))

(defn f-deriv-2 [x] ; second derivative
  ($= 6 * x))

(def i-vals [0.5 1.5])

(def step-sizes [0.00001 0.0001 0.001 0.01 0.1 0.2 0.3])

(def meths [:central :backward :forward])

;; Functions to calculate aggregate data for 2.12 ;;
(defn calc-error
  "Calculates percentage error"
  [actual expected]
  ($= 100 * (actual - expected) / expected))

(defn run
  "Runs through all scenarios for problem 2.12.
   Splits up data by x-value and method used,
   aggregating approximations for all step sizes."
  []
  (doseq [power [1 2]]
    (let [title-power (case power     ; set power component of title
                        1 "1st-deriv"
                        2 "2nd-deriv")
          deriv-func  (case power     ; set derivative function based on power
                        1 f-deriv
                        2 f-deriv-2)]
    (doseq [x i-vals]   ; iterate through x values (0.5 1.5)
      (doseq [m meths]  ; iterate through methods (central, backward, forward)
        (let [title (str x "-" (name m) "-" title-power)
              data  (for [step step-sizes]                                ; iterate through x-steps
                      (let [y       (fd-approx power 2 m f x step)        ; compute at x with given method, step-size
                            error   (abs (calc-error y (deriv-func x)))]  ; compute error
                         {:step-size  step
                         :value       y
                         :error       error}))]
          (save (to-dataset data) (str "/users/jadis/documents/cfd/chapter-2/" title ".csv"))))))))
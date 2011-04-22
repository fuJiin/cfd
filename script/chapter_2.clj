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
  (doseq [x i-vals]   ; iterate through x values (0.5 1.5)
    (doseq [m meths]  ; iterate through methods (central, backward, forward)
      (let [title (str x "-" (name m))
            data  (for [step step-sizes]                              ; iterate through x-steps
                    (let [y       (fd-approx 2 2 m f x step)          ; compute at x with given method, step-size
                          error   (abs (calc-error y (f-deriv-2 x)))] ; compute error
                       {:step-size  step
                       :value       y
                       :error       error}))]
        (save (to-dataset data) (str "/users/jadis/desktop/" title ".csv"))))))
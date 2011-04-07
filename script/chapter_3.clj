(ns cfd.script.chapter-3
  (:use [incanter.core])
  (:use [cfd.odppde]))

(defn analytic-solve
  "Analytical solution for 3.1"
  [ts ti l x t alpha n]
  (let [f (fn [m]
            ($= (Math/E ** (0 - (((m * Math/PI / l) ** 2) * alpha * t)))
                *
                (1 - (-1 ** m))
                *
                (sin ($= m * Math/PI * x / l))
                /
                (m * Math/PI)))
        series (reduce +
                (take n
                  (pmap f (iterate inc 1))))]
    ($= ts + 2 * (ti - ts) * series)))
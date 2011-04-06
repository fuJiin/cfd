(ns cfd.script.chapter-3
  (:use [incanter.core])
  (:use [cfd.odppde]))
  
(defn analytic
  [ti ts l x n alpha t]
  (let [f (fn [m]
            ($= (Math/E ** (0 - (((m * Math/PI / l) ** 2) * alpha * t)))
                *
                (1 - (-1 ** m))
                *
                (sin (m * Math/PI * x / l))
                /
                (m * Math/PI)))
        series (take n
                (map f (iterate inc 1)))]
    ($= ti + 2 * (ti - ts) * series)))
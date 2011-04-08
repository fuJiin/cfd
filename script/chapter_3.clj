(ns cfd.script.chapter-3
  (:use [incanter core charts])
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

;; Set problem conditions ;;
(def thickness 1)     ; ft        ; wall thickness
(def temp-init 100)   ; deg F     ; initial temperature
(def temp-surf 300)   ; deg F     ; surface temperature
(def alpha 0.1)       ; ft^2 / hr ; diffusivity
(def t-limit 0.5)     ; hrs       ; time limit
(def display-int 0.1) ; hr        ; display-interval

;; Graph helpers ;;
(defn init-grid
  "Initialize grid based on x-steps, surface temperature,
   and initial temperature."
   [x-step]
   (let [x-steps ($= thickness / x-step + 1)] ; add 1 to account for starting pt
     (vec
      (for [i (range 0 x-steps)]
        (if (or (= i 0) (= i (dec x-steps)))
            temp-surf
            temp-init)))))

(defn simulate
  "Simulate scenario given an approximation technique f,
   x-step, and t-step"
  [f x-step t-step]
  (let [t-steps   ($= t-limit / t-step)
        x-steps   ($= thickness / x-step + 1) ; add 1 to account for starting pt
        grid      (init-grid x-step)]
    (if (= f dufort-frankel)
        (dufort-frankel-run t-steps grid x-step t-step alpha)
        (regular-run f t-steps grid x-step t-step alpha))))

(defn get-disp-data
  "Fetch display data from simulation, using the display interval.
   This allows t-step and display-t-interval to be different."
  [f x-step t-step]
  (let [data        (simulate f x-step t-step)
        t-steps     ($= t-limit / t-step)
        disp-steps  ($= t-steps * display-int / t-limit)
        t-indices   (map #($= % * disp-steps)
                          (take ($= t-steps / disp-steps + 1) (iterate inc 0)))]
    (vec
      (for [i t-indices]
        (nth data i)))))

(defn make-dataset
  "Transform display data into an incanter dataset."
  [f x-step t-step]
  (let [data        (get-disp-data f x-step t-step)
        t-steps     ($= t-limit / t-step)
        t-column    (map #(* % t-step)
                          (take (inc t-steps) (iterate inc 0)))
        disp-steps  ($= t-steps * display-int / t-limit)
        x-steps     ($= thickness / x-step)
        x-column    (map #(* % x-step)
                          (take (inc x-steps) (iterate inc 0)))]
    (loop [row-i 0
           d-set []]
      (if (= row-i (count data))
          (to-dataset d-set)
          (let [row     (nth data row-i)
                mapping (for [i (range 0 (count row))]
                          {"x" (nth x-column i)
                           "T" (nth row i)
                           "t" (nth t-column row-i)})]
            (recur (inc row-i)
                   (into d-set mapping)))))))

(defn graph
  "Graph simulation dataset"
  [f x-step t-step]
  (let [d-set (make-dataset f x-step t-step)]
    (with-data d-set
      (line-chart :x :T
        :title    (str ((meta f) :name)) ; name from function metadata
        :group-by :t
        :legend   true))))
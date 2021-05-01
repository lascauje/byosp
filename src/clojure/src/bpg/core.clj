"Implementation of Convex-Concave Backtracking for Inertial
Bregman Proximal Gradient Algorithms in Non-Convex Optimization."

(ns bpg.core
  (:require [uncomplicate.neanderthal
             [core :refer [nrm1 nrm2 mrows ncols dim axpy trans mm xpy scal]]
             [native :refer [dge]]
             [linalg :refer [ev!]]
             [math :refer [pow]]
             [real :refer [entry]]
             [vect-math :refer [abs fmax copy-sign]]]))

(defn mx->scal
  "Convert a matrix 'a' of dim 1x1 to a scalar."
  [a]
  (entry a 0 0))

(defn scal->mx
  "Convert a scalar 's' to a matrix 'a' of dim 'a'."
  [s a]
  (dge (mrows a) (ncols a) (repeat (dim a) s)))

(defn mxpy
  "Substract vectors or matrices."
  [a b]
  (axpy -1 b a))

(defn roots
  "Return the roots of a polynomial with coefficients given in p.

  p[0] * x**n + p[1] * x**(n-1) + ... + p[n-1]*x + p[n]

  The algorithm relies on computing the eigenvalues of the companion matrix.
  See numpy.roots

  NB: p must be only of degree 3."
  [p]
  (ev! (dge 3 3 [0 0 (/ (p 3) (- (p 0)))
                 1 0 (/ (p 2) (- (p 0)))
                 0 1 (/ (p 1) (- (p 0)))])))

(defn f
  "See p.20 (6.5)."
  [lambda x]
  (* lambda (nrm1 x)))

(defn g
  "See p.19 (6.5)."
  [a x b]
  (letfn [(g-helper
            [a x b]
            (* 0.25 (pow (- (mx->scal (mm (trans x) (mm a x))) b) 2)))]
    (apply + (for [i (range (count a))]
               (g-helper (a i) x (b i))))))
(defn d-g
  "See p.19 (6.5)."
  [a x b]
  (letfn [(d-g-helper
            [a x b]
            ;; It should be "b**2, but the official implementation uses "b"
            (scal (- (mx->scal (mm (trans x) (mm a x))) b)
                  (mm a x)))]
    (apply xpy (for [i (range (count a))]
                 (d-g-helper (a i) x (b i))))))

(defn h
  "See p.20 (6.6)."
  [x]
  (+ (* 0.25 (pow (nrm2 x) 4))
     (* 0.5 (pow (nrm2 x) 2))))

(defn d-h
  "See p.20 (6.6)."
  [x]
  (xpy (scal (mx->scal (mm (trans x) x)) x) x))

(defn obj-func
  "See p.19 (6.5)."
  [a x b lambda]
  (+ (f lambda x) (g a x b)))

(defn breg-dist
  "See p.4 definition 2.1."
  [x y]
  (let [dist (- (h x) (h y) (mx->scal (mm (trans (mxpy x y)) (d-h y))))]
    (if (< 1e-15 dist)
      dist
      0)))

(defn find-y
  "See p.7 (3.4)."
  [gamma prev-x x]
  (xpy x (scal gamma (mxpy x prev-x))))

(defn find-gamma
  "See p.7 (3.5)."
  ([delta epsilon gamma prev-x x lower upper prev-tau]
   (let [y (find-y gamma prev-x x)]
     (if (< (/ (* (- delta epsilon) (breg-dist prev-x x))
               (inc (* lower prev-tau)))
            (breg-dist x y))
       (recur delta epsilon (* gamma 0.9) prev-x x lower upper prev-tau)
       [y gamma])))

  ([delta epsilon prev-x x lower upper]
   (find-gamma delta epsilon 1 prev-x x lower upper (/ 1 upper))))

(defn find-lower
  "See p.7 (3.5)."
  [delta epsilon a b prev-x x lower upper]
  (let [[y gamma] (find-gamma delta epsilon prev-x x lower upper)]
    (if (< 1e-7 (- (+ (g a y b) (mx->scal (mm (trans (d-g a y b)) (mxpy x y))))
                   (* lower (breg-dist x y))
                   (g a x b)))
      ;; Scaling "lower L" parameter with 2
      (recur delta epsilon a b prev-x x (* lower 2) upper)
      [lower y gamma])))

(defn soft-thresh
  "See p.21 (6.10)."
  [theta y]
  (copy-sign (fmax (mxpy (abs y) (scal->mx theta y)) (scal->mx 0 y)) y))

(defn t-star
  "See p.21 (6.10) next."
  [coeff]
  (entry (roots [(pow (nrm2 coeff) 2) 0 1 -1]) 0 0))

(defn find-x
  "See both p.7 (3.7) and p.21 (6.10) next."
  [lambda a b y upper]
  (let [tau (/ 1 upper)
        s (soft-thresh (* lambda tau) (mxpy (d-h y) (scal tau (d-g a y b))))]
    (scal (t-star s) s)))

(defn find-upper
  "See p.7 (3.5)."
  [lambda a b y lower upper]
  (let [x (find-x lambda a b y upper)]
    (if (< (- (+ (g a y b)
                 (mx->scal (mm (trans (d-g a y b)) (mxpy x y)))
                 (* upper (breg-dist x y)))
              (g a x b))
           -1e-7)
      ;; Scaling "upper L" parameter with 2
      (recur lambda a b y lower (* upper 2))
      [upper x])))

(defn find-upper-official
  "See p.7 (3.5).

    This is a wrapper to the raw find-upper function.
    In the raw find-upper function the paper has been implemented
    slightly differently, but in order to check the results
    this function (find-upper-official) is used instead."
  [lambda a b y lower upper]
  (let [x (find-x lambda a b y upper)]
    (if (< (- (+ (g a y b)
                 (mx->scal (mm (trans (d-g a y b)) (mxpy x y)))
                 (* upper (breg-dist x y)))
              (g a x b))
           -1e-7)
      (let [[upper _] (find-upper lambda a b y lower upper)
            next-upper (* upper 2)
            next-x (find-x lambda a b y next-upper)]
        [next-upper next-x])
      [upper x])))

(defn lyapunov
  "See p.11 (5.3) (simpler)."
  [lambda a b prev-x x upper]
  (let [prev-tau (/ 1 upper)]
    (+ (* prev-tau (obj-func a x b lambda)) (breg-dist x prev-x))))

(defn coca-bpg
  "See p.7 Convex-Concave Inertial BPG algorithm."
  [lambda delta epsilon a b prev-x x prev-lower prev-upper]
  (let [[lower y gamma] (find-lower delta epsilon a b prev-x x
                                    prev-lower prev-upper)
        [upper next-x] (find-upper-official lambda a b y lower prev-upper)
        res (lyapunov lambda a b x next-x upper)]
    [gamma x next-x y lower upper res]))

(defn iter-coca-bpg
  "Wrapper to iterate on coca-bpg function."
  [lambda delta epsilon a b prev-x x prev-lower prev-upper nb-iter]
  (loop [i 0
         [gamma x next-x y lower upper res] (coca-bpg lambda delta epsilon
                                                      a b prev-x x
                                                      prev-lower prev-upper)
         res-all []]
    (if (< i (dec nb-iter))
      (recur (inc i)
             (coca-bpg lambda delta epsilon a b x next-x lower upper)
             (conj res-all res))
      [lambda delta epsilon gamma
       a b x next-x y lower upper (conj res-all res)])))

(defn lcg
  "A simple linear congruential generator inspired from the glibc.

  Values are in range ]0,1[.
  See glibc/stdlib/random_r.c:__random_r:341"
  [nb-time seed]
  (letfn [(iterator [a b]
            (fn [x] (mod (+ (* a x) b) (bit-shift-left 1 31))))]
    (let [glibc (drop 1 (iterate (iterator 1103515245 12345) seed))]
      (mapv #(double (/ % (bit-shift-left 1 31))) (take nb-time glibc)))))

(defn init-a
  "Init of a (list of matrices) using a linear congruential generator.

    See p.19 (6.3)."
  [nb-time row]
  (letfn  [(single-matrix [seed]
             (let [m (dge row 1 (lcg row seed))]
               (mm m (trans m))))]
    (mapv single-matrix (range nb-time))))

(defn init-b
  "Init of b (list of float) using a linear congruential generator.

    See p.19 (6.3)."
  [nb-time]
  (lcg nb-time 0))

(defn bench
  "Benchmark the BPG algorithm."
  ([dim nb-time nb-iter]
   (let [lambda 0.1
         delta 0.15
         epsilon 0.00001
         a (init-a nb-time dim)
         b (init-b nb-time)
         prev-x (dge dim 1 (repeat dim 1))
         x (dge dim 1 (repeat dim 1))
         prev-lower 0.001
         prev-upper 10]
     (time (last (last (iter-coca-bpg lambda
                                      delta
                                      epsilon
                                      a
                                      b
                                      prev-x
                                      x
                                      prev-lower
                                      prev-upper
                                      nb-iter))))))
  ([dim]
   (bench dim 100 1000)))

(defn -main
  "I don't do a whole lot ... yet."
  []
  (println "Hello BPG!"))

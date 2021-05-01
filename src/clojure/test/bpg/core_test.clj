"Unit tests for Convex-Concave Backtracking for Inertial
Bregman Proximal Gradient Algorithms in Non-Convex Optimization."

(ns bpg.core-test
  (:require [clojure.test :refer :all]
            [bpg.core :refer :all]
            [uncomplicate.neanderthal
             [native :refer :all]
             [core :refer [sum]]
             [math :refer [f=]]]))

(deftest roots-test
  (testing "Sanity check."
    (is (f= (sum (roots [10, 0, 1, -1]))
            -5.551115123125783E-17))))

(deftest f-test
  (testing "Sanity check."
    (is (f= (f 0.1 (dge 2 1 [1 1]))
            0.2))))

(deftest g-test
  (testing "Sanity check."
    (is (f= (g [(dge 2 2 [1 2 3 4]) (dge 2 2 [1 2 3 4])]
               (dge 2 1 [1 2])
               [0.60276 0.64589])
            347.838348377425))))

(deftest d-g-test
  (testing "Sanity check."
    (is (f= (sum (d-g [(dge 2 2 [1 2 3 4]) (dge 2 2 [1 2 3 4])]
                      (dge 2 1 [1 2])
                      [0.60276 0.64589]))
            896.77295))))

(deftest h-test
  (testing "Sanity check."
    (is (f= (h (dge 2 1 [1 2]))
            8.75))))

(deftest d-h-test
  (testing "Sanity check."
    (is (f= (sum (d-h (dge 2 1 [1 2])))
            18.0))))

(deftest obj-func-test
  (testing "Sanity check."
    (is (f= (obj-func [(dge 2 2 [1 2 3 4]) (dge 2 2 [1 2 3 4])]
                      (dge 2 1 [1 2])
                      [0.60276 0.64589]
                      0.1)
            348.138348377425))))

(deftest breg-dist-test
  (testing "Sanity check."
    (is (f= (breg-dist (dge 2 1 [3 4]) (dge 2 1 [1 2]))
            124.0))))

(deftest find-y-test
  (testing "Sanity check."
    (is (f= (sum (find-y 0.9 (dge 2 1 [3 4]) (dge 2 1 [10 20])))
            50.7))))

(deftest find-gamma-test
  (testing "Sanity check."
    (let [[y gamma] (find-gamma 0.15 0.00001
                                (dge 2 1 [3 4]) (dge 2 1 [10 20])
                                0.001 10)]
      (is (and (f= (sum y) 35.846291405156705)
               (f= gamma 0.2541865828329001))))))

(deftest soft-thresh-test
  (testing "Sanity check."
    (is (f= (sum (soft-thresh 0.01 (dge 2 1 [-1.92486 0.15026])))
            -1.7746))))

(deftest t-star-test
  (testing "Sanity check."
    (is (f= (t-star (dge 2 1 [-1.92486 0.15026]))
            0.5088523346606675))))

(deftest find-x-test
  (testing "Sanity check."
    (is (f= (sum (find-x 0.1
                         [(dge 2 2 [1 2 3 4]) (dge 2 2 [1 2 3 4])]
                         [0.60276 0.64589]
                         (dge 2 1 [-1 0])
                         10))
            -0.9054999300267382))))

(deftest lyapunov-test
  (testing "Sanity check."
    (is (f= (lyapunov 0.1
                      [(dge 2 2 [1 2 3 4]) (dge 2 2 [1 2 3 4])]
                      [0.60276 0.64589]
                      (dge 2 1 [3 4])
                      (dge 2 1 [1 2])
                      10)
            238.8138348377425))))

(deftest find-lower-test
  (testing "Sanity check."
    (let [[lower y gamma] (find-lower 0.15
                                      0.00001
                                      [(dge 2 2 [0.99531367 0.9877355
                                                 0.58327082 0.66336516])
                                       (dge 2 2 [0.0981552 0.61172762
                                                 0.52265855 0.0371408])]
                                      [0.6027633760716439 0.6458941130666561]
                                      (dge 2 1 [0.4719454622448402
                                                0.4220473240632225])
                                      (dge 2 1 [0.4693035792040524
                                                0.42080845789332616])
                                      0.001
                                      10)]
      (is (and (f= lower 0.512)
               (f= (sum y) 0.8887589035161779)
               (f=  gamma 0.34867844010000015))))))

(deftest find-upper-test
  (testing "Sanity check."
    (let [[upper x] (find-upper 0.1
                                [(dge 2 2 [0.99531367 0.9877355
                                           0.58327082 0.66336516])
                                 (dge 2 2 [0.0981552 0.61172762
                                           0.52265855 0.0371408])]
                                [0.6027633760716439 0.6458941130666561]
                                (dge 2 1 [0.4588713468854095
                                          0.4191807008769411])
                                32.768
                                10)]
      (is (and (f= upper 80.0)
               (f= (sum x) 0.8779832567675403))))))

(deftest find-upper-official-test
  (testing "Sanity check."
    (let [[upper x] (find-upper-official 0.1
                                         [(dge 2 2 [0.99531367 0.9877355
                                                    0.58327082 0.66336516])
                                          (dge 2 2 [0.0981552 0.61172762
                                                    0.52265855 0.0371408])]
                                         [0.6027633760716439
                                          0.6458941130666561]
                                         (dge 2 1 [0.4588713468854095
                                                   0.4191807008769411])
                                         32.768
                                         10)]
      (is (and (f= upper 160.0)
               (f= (sum x) 0.8780176535863605))))))

(deftest lcg-test
  (testing "Sanity check."
    (is (f= (apply + (lcg 10 0))
            3.977038038428873))))

(deftest iter-coca-bpg-test-10
  (testing "Sanity check. Lyapunov values must be decreasing."
    (let [dim 10
          nb-time 100
          nb-iter 1000
          lambda 0.1
          delta 0.15
          epsilon 0.00001
          a (init-a nb-time dim)
          b (init-b nb-time)
          prev-x (dge dim 1 (repeat dim 1))
          x (dge dim 1 (repeat dim 1))
          prev-lower 0.001
          prev-upper 10
          res (iter-coca-bpg
               lambda delta epsilon a b prev-x x prev-lower prev-upper
               nb-iter)]
      (is (and (apply >= (last res))
               (f= (last (last res)) 9.734142808925895E-4))))))

(deftest iter-coca-bpg-test-128
  (testing "Sanity check. Lyapunov values must be decreasing."
    (let [dim 128
          nb-time 100
          nb-iter 1000
          lambda 0.1
          delta 0.15
          epsilon 0.00001
          a (init-a nb-time dim)
          b (init-b nb-time)
          prev-x (dge dim 1 (repeat dim 1))
          x (dge dim 1 (repeat dim 1))
          prev-lower 0.001
          prev-upper 10
          res (iter-coca-bpg
               lambda delta epsilon a b prev-x x prev-lower prev-upper
               nb-iter)]
      (is (and (apply >= (last res))
               (f= (last (last res)) 8.69888363198627E-5))))))

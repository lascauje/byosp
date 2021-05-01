(defproject bpg "0.1.0"
  :description "Convex-Concave Backtracking for Bregman Proximal Gradient"
  :url "https://paperswithcode.com/paper/convex-concave-backtracking-for-inertial"
  :license {:name "GPL-3.0"
            :url "https://www.gnu.org/licenses/gpl-3.0.en.html"}
  :dependencies [[org.clojure/clojure "1.10.1"]
                 [uncomplicate/neanderthal "0.41.0"]
                 [org.bytedeco/mkl-platform-redist "2020.3-1.5.4"]
                 [clj-kondo "2021.03.31"]
                 [codox-theme-rdash "0.1.2"]]
  :plugins [[jonase/eastwood "0.4.0"]
            [lein-kibit "0.1.8"]
            [lein-cljfmt "0.7.0"]
            [lein-codox "0.10.7"]]
  :codox {:metadata {:doc/format :markdown}
          :themes [:rdash]}
  :aliases {"kondo" ["run" "-m" "clj-kondo.main"]
            "ci" ["do"
                  ["eastwood"]
                  ["kondo" "--lint" "src"]
                  ["cljfmt" "check"]
                  ["kibit"]
                  ["test"]
                  ["codox"]]}
  :jvm-opts ^:replace ["--add-opens=java.base/jdk.internal.ref=ALL-UNNAMED"]
  :main ^:skip-aot bpg.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all
                       :jvm-opts ["-Dclojure.compiler.direct-linking=true"]}})

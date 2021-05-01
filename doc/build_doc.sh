#!/bin/bash -eu

rm -rf _build/
sphinx-build . _build/
cd ../src/clojure/ && lein codox && cd -
cp -R ../src/clojure/target/default/doc/ _build/clj
cd ../src/python/ && nox -rs docs && cd -
cp -R ../src/python/docs/_build/ _build/py
rm -rf ../bin/
mkdir ../bin/ && mv _build/ ../bin/doc

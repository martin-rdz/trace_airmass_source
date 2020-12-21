#!/bin/bash

# generate the docs and open them in firefox
cd docs
make html
firefox /home/larda/larda-doc/html/index.html

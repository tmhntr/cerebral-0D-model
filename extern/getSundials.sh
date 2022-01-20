#!/bin/bash

curl https://github.com/LLNL/sundials/releases/download/v2.7.0/sundials-2.7.0.tar.gz -o sundials-2.7.0.tar.gz -L
tar xzf sundials-2.7.0.tar.gz 
mkdir builddir
cd builddir

cmake ../sundials-2.7.0
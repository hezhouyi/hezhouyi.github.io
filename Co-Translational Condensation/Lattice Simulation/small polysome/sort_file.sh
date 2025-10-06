#!/bin/bash

mkdir -p log
mv *.err log
mv *.out log
mv *.slurm log

mkdir -p RawData
mv phase* RawData

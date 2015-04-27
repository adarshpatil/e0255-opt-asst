#!/bin/bash
cat ref | sort -k1,1n -k2,2n -o ref1 &
cat opt | sort -k1,1n -k2,2n -o opt1 &
wait
head -n 500 ref1 > ref500
head -n 500 opt1 > opt500

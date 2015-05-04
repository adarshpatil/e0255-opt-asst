#!/bin/bash
if [  $# -ne 1 ]
then
echo "Usage: $0 <log-file-name>"
fi
split -l 50 $1
echo "Average OpenCV Time"
cat xaa | awk -F" " '{print $4}' | paste -sd+ | bc | awk '{print $1/50}'
echo "Average Reference Time"
cat xab | awk -F" " '{print $4}' | paste -sd+ | bc | awk '{print $1/50}'
echo "Average Optimized Time"
cat xac | awk -F" " '{print $4}' | paste -sd+ | bc | awk '{print $1/50}'

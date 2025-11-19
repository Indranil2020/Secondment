#!/bin/bash
# GPU4PySCF Runner Script
# Fixes CUDA library path issue

# Set CUDA 12.9 library path
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH

# Run the GPU-accelerated script
python3 tdm_calc_accurate_GPU.py "$@"

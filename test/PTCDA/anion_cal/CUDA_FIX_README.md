# CUDA Library Path Fix for GPU4PySCF

## 🔍 The Problem

```
ImportError: /usr/local/cuda/targets/x86_64-linux/lib/libcublas.so.12: 
undefined symbol: cublasLtGetEnvironmentMode, version libcublasLt.so.12
```

**Root Cause:**
- CUDA Toolkit: 11.5 (from `nvcc`)
- CUDA Driver: 12.9 (from `nvidia-smi`)
- GPU4PySCF: Installed for CUDA 12.x
- System default `/usr/local/cuda` points to CUDA 11.5
- GPU4PySCF needs CUDA 12.9 libraries

## ✅ Quick Fix (Temporary)

Run with the wrapper script:
```bash
./run_gpu.sh
```

Or manually set the path:
```bash
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH
python3 tdm_calc_accurate_GPU.py
```

## ✅ Permanent Fix (Recommended)

### Option 1: Add to ~/.bashrc (For Your User Only)

```bash
echo 'export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

### Option 2: Update CUDA Symlink (System-Wide, Requires Root)

```bash
sudo rm /usr/local/cuda
sudo ln -s /usr/local/cuda-12.9 /usr/local/cuda
```

### Option 3: Use the Wrapper Script

Always run:
```bash
./run_gpu.sh
```

## 🧪 Verify the Fix

Test GPU4PySCF import:
```bash
python3 -c "from gpu4pyscf import dft; print('✓ GPU4PySCF OK')"
```

Should output:
```
✓ GPU4PySCF OK
```

## 📊 System Information

Your system has:
- **CUDA Toolkit 11.5:** `/usr/local/cuda` → `/etc/alternatives/cuda`
- **CUDA Runtime 12.9:** `/usr/local/cuda-12.9/` (needed for GPU4PySCF)
- **NVIDIA Driver:** 575.51.03 (supports CUDA 12.9)
- **GPU4PySCF:** 1.4.3 (CUDA 12.x version)

## 🎯 Recommended Action

**For permanent fix, add to ~/.bashrc:**

```bash
# Add CUDA 12.9 libraries for GPU4PySCF
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH
```

Then:
```bash
source ~/.bashrc
python3 tdm_calc_accurate_GPU.py  # Should work now!
```

## 🚀 Usage

### With Wrapper Script (No Setup Needed)
```bash
./run_gpu.sh
```

### After Permanent Fix
```bash
python3 tdm_calc_accurate_GPU.py
```

### CPU Version (No CUDA Issues)
```bash
python3 tdm_calc_accurate.py
```

## ⚠️ Why This Happened

1. Your system has **two CUDA versions** installed
2. Default `/usr/local/cuda` points to **CUDA 11.5**
3. GPU4PySCF needs **CUDA 12.x** libraries
4. Without `LD_LIBRARY_PATH`, Python loads wrong CUDA version
5. Library version mismatch causes the symbol error

## 🔧 Alternative: Reinstall for CUDA 11.x

If you prefer to use CUDA 11.5:
```bash
pip uninstall gpu4pyscf-cuda12x
pip install gpu4pyscf-cuda11x
```

But CUDA 12.9 is newer and recommended!

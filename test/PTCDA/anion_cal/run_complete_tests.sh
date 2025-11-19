#!/bin/bash
# Complete test suite: H2O with charges 0, +1, -1 for CPU and GPU

echo "=========================================="
echo "COMPLETE TEST SUITE: H2O (CPU vs GPU)"
echo "=========================================="
echo ""

# Test 1: CPU H2O Neutral (charge=0)
echo "Test 1/6: CPU H2O Neutral (charge=0)..."
sed -i 's/^CHARGE = .*/CHARGE = 0/' tdm_calc_accurate.py
python3 tdm_calc_accurate.py > test_final_cpu_charge0.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Test 1 PASSED"
else
    echo "❌ Test 1 FAILED"
fi
echo ""

# Test 2: CPU H2O Cation (charge=+1)
echo "Test 2/6: CPU H2O Cation (charge=+1)..."
sed -i 's/^CHARGE = .*/CHARGE = 1/' tdm_calc_accurate.py
python3 tdm_calc_accurate.py > test_final_cpu_charge+1.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Test 2 PASSED"
else
    echo "❌ Test 2 FAILED"
fi
echo ""

# Test 3: CPU H2O Anion (charge=-1)
echo "Test 3/6: CPU H2O Anion (charge=-1)..."
sed -i 's/^CHARGE = .*/CHARGE = -1/' tdm_calc_accurate.py
python3 tdm_calc_accurate.py > test_final_cpu_charge-1.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Test 3 PASSED"
else
    echo "❌ Test 3 FAILED"
fi
echo ""

# Test 4: GPU H2O Neutral (charge=0)
echo "Test 4/6: GPU H2O Neutral (charge=0)..."
sed -i 's/^CHARGE = .*/CHARGE = 0/' tdm_calc_accurate_GPU.py
./run_gpu.sh > test_final_gpu_charge0.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Test 4 PASSED"
else
    echo "❌ Test 4 FAILED"
fi
echo ""

# Test 5: GPU H2O Cation (charge=+1)
echo "Test 5/6: GPU H2O Cation (charge=+1)..."
sed -i 's/^CHARGE = .*/CHARGE = 1/' tdm_calc_accurate_GPU.py
./run_gpu.sh > test_final_gpu_charge+1.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Test 5 PASSED"
else
    echo "❌ Test 5 FAILED"
fi
echo ""

# Test 6: GPU H2O Anion (charge=-1)
echo "Test 6/6: GPU H2O Anion (charge=-1)..."
sed -i 's/^CHARGE = .*/CHARGE = -1/' tdm_calc_accurate_GPU.py
./run_gpu.sh > test_final_gpu_charge-1.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Test 6 PASSED"
else
    echo "❌ Test 6 FAILED"
fi
echo ""

echo "=========================================="
echo "EXTRACTING RESULTS..."
echo "=========================================="
echo ""

# Extract and compare results
python3 << 'EOF'
import re

tests = [
    ("CPU charge=0", "test_final_cpu_charge0.log"),
    ("CPU charge=+1", "test_final_cpu_charge+1.log"),
    ("CPU charge=-1", "test_final_cpu_charge-1.log"),
    ("GPU charge=0", "test_final_gpu_charge0.log"),
    ("GPU charge=+1", "test_final_gpu_charge+1.log"),
    ("GPU charge=-1", "test_final_gpu_charge-1.log"),
]

print(f"{'Test':<20} {'Spin':<6} {'Method':<6} {'Energy (a.u.)':<15} {'State 1 (eV)':<12}")
print("="*70)

for name, logfile in tests:
    try:
        with open(logfile, 'r') as f:
            content = f.read()
            
        spin_match = re.search(r'Spin multiplicity \(2S\+1\): (\d+)', content)
        method_match = re.search(r'DFT method: (\w+)', content)
        energy_match = re.search(r'Ground state energy: ([-\d.]+)', content)
        state1_match = re.search(r'State 1: [\d.]+ a\.u\. = ([\d.]+) eV', content)
        
        spin = spin_match.group(1) if spin_match else "N/A"
        method = method_match.group(1) if method_match else "N/A"
        energy = energy_match.group(1) if energy_match else "N/A"
        state1 = state1_match.group(1) if state1_match else "N/A"
        
        print(f"{name:<20} {spin:<6} {method:<6} {energy:<15} {state1:<12}")
    except:
        print(f"{name:<20} {'ERR':<6} {'ERR':<6} {'ERROR':<15} {'ERROR':<12}")

print("="*70)
EOF

echo ""
echo "=========================================="
echo "TEST SUITE COMPLETE!"
echo "=========================================="

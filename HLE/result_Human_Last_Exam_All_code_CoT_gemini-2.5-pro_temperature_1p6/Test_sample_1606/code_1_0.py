import math

# Given parameters
# The linear convolution is between a length-1200 sequence and a length-90 sequence.
# For fast convolution algorithms like overlap-add/save, we segment the longer sequence.
L = 1200  # Length of the long sequence (input signal x[n])
M = 90    # Length of the short sequence (filter h[n])
N = 128   # N-point DFT/IDFT size

print("This script calculates the number of DFT and IDFT operations for Overlap-Add vs. Overlap-Save methods.\n")

# --- Overlap-Add Method Calculation ---
print("### Overlap-Add Method ###\n")

# Step 1: Determine the length of the data blocks.
# The constraint is N >= L_block + M - 1. We choose the maximum possible L_block.
L_block_add = N - M + 1

# Step 2: Determine the number of blocks needed.
# The number of blocks K is ceil(L / L_block).
K_add = math.ceil(L / L_block_add)

# Step 3: Calculate the total number of DFT and IDFT operations.
num_dfts_add = K_add + 1  # K for input blocks + 1 for the filter
num_idfts_add = K_add       # K for the output blocks
total_ops_add = num_dfts_add + num_idfts_add

# Print the step-by-step calculation
print(f"1. Determine data block length (L_block):")
print(f"   L_block = N - M + 1 = {N} - {M} + 1 = {L_block_add}")

print(f"\n2. Determine the number of blocks (K):")
print(f"   K = ceil(L / L_block) = ceil({L} / {L_block_add}) = {K_add}")

print(f"\n3. Calculate total operations:")
print(f"   Number of DFTs = K (for input blocks) + 1 (for filter) = {K_add} + 1 = {num_dfts_add}")
print(f"   Number of IDFTs = K (for output blocks) = {K_add}")
print(f"   Total (DFT + IDFT) operations = {num_dfts_add} + {num_idfts_add} = {total_ops_add}")

print("\n" + "="*40 + "\n")

# --- Overlap-Save Method Calculation ---
print("### Overlap-Save Method ###\n")

# Step 1: Determine the number of new data points per block.
# Each block produces N - M + 1 valid output samples. This corresponds to the number of new input samples.
L_useful_save = N - M + 1

# Step 2: Determine the number of blocks needed.
# The number of blocks K is ceil(L / L_useful).
K_save = math.ceil(L / L_useful_save)

# Step 3: Calculate the total number of DFT and IDFT operations.
num_dfts_save = K_save + 1  # K for input blocks + 1 for the filter
num_idfts_save = K_save       # K for the output blocks
total_ops_save = num_dfts_save + num_idfts_save

# Print the step-by-step calculation
print(f"1. Determine new data points per block:")
print(f"   New points = N - M + 1 = {N} - {M} + 1 = {L_useful_save}")

print(f"\n2. Determine the number of blocks (K):")
print(f"   K = ceil(L / New points) = ceil({L} / {L_useful_save}) = {K_save}")

print(f"\n3. Calculate total operations:")
print(f"   Number of DFTs = K (for input blocks) + 1 (for filter) = {K_save} + 1 = {num_dfts_save}")
print(f"   Number of IDFTs = K (for output blocks) = {K_save}")
print(f"   Total (DFT + IDFT) operations = {num_dfts_save} + {num_idfts_save} = {total_ops_save}")
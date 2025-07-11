import math

# Given parameters for the convolution problem
L = 90    # Length of the first sequence x[n] (the filter)
M = 1200  # Length of the second sequence h[n] (the signal)
N = 128   # DFT/IDFT size

print("--- Calculating the number of DFT/IDFT operations ---")
print(f"Parameters: L={L}, M={M}, N={N}")
print("-" * 50)

# --- Overlap-Add Method ---
print("Overlap-Add Method:")

# Step 1: Determine the optimal block length for the signal h[n].
# To avoid time-domain aliasing in linear convolution, the DFT size N must be
# greater than or equal to the length of the convolved block (L + L_b - 1).
# We choose the maximum possible block length L_b to minimize the number of blocks.
# L_b = N - L + 1
L_b_add = N - L + 1
print(f"1. Block length (L_b) = N - L + 1 = {N} - {L} + 1 = {L_b_add}")

# Step 2: Calculate the number of blocks (K) needed to process the entire signal.
# This is the ceiling of the total signal length divided by the block length.
K_add = math.ceil(M / L_b_add)
print(f"2. Number of blocks (K) = ceil(M / L_b) = ceil({M} / {L_b_add}) = {K_add}")

# Step 3: Calculate the total number of (DFT+IDFT) operations.
# This includes 1 DFT for the filter x[n], K DFTs for the signal blocks,
# and K IDFTs for the output blocks.
dfts_add = K_add + 1
idfts_add = K_add
total_add = dfts_add + idfts_add
print(f"3. Total operations = (K+1 DFTs) + (K IDFTs) = ({K_add}+1) + {K_add} = {total_add}")
print("-" * 50)


# --- Overlap-Save Method ---
print("Overlap-Save Method:")

# Step 1: Determine the number of valid output samples per block.
# Each N-point circular convolution produces N - (L-1) valid samples.
L_data_save = N - L + 1
print(f"1. Valid samples per block = N - L + 1 = {N} - {L} + 1 = {L_data_save}")

# Step 2: Calculate the number of blocks (K) needed.
# This is the ceiling of the total signal length divided by the valid samples per block.
K_save = math.ceil(M / L_data_save)
print(f"2. Number of blocks (K) = ceil(M / L_data) = ceil({M} / {L_data_save}) = {K_save}")

# Step 3: Calculate the total number of (DFT+IDFT) operations.
# This includes 1 DFT for the filter x[n], K DFTs for the signal blocks,
# and K IDFTs for the output blocks.
dfts_save = K_save + 1
idfts_save = K_save
total_save = dfts_save + idfts_save
print(f"3. Total operations = (K+1 DFTs) + (K IDFTs) = ({K_save}+1) + {K_save} = {total_save}")
print("-" * 50)

print(f"Final Answer: The number of (DFT+IDFT) operations is {total_add} for the overlap-add method and {total_save} for the overlap-save method.")
<<<63 for overlap add and 63 for overlap save>>>
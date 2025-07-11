import math

# Step 1: Define Parameters
# The lengths of the two sequences and the DFT size are given.
L_x = 1200  # Length of the long sequence
L_h = 90   # Length of the short sequence (filter)
N = 128    # DFT/IDFT size

print("--- Problem Parameters ---")
print(f"Length of long sequence (L_x): {L_x}")
print(f"Length of short sequence (L_h): {L_h}")
print(f"DFT/IDFT size (N): {N}")
print("-" * 30)

# Step 2: Calculate operations for Overlap-Add method
print("--- Overlap-Add Method Calculation ---")
# To prevent time-domain aliasing in circular convolution, the DFT size N must be
# greater than or equal to the linear convolution length of one block and the filter.
# N >= L + L_h - 1, where L is the block size.
# To maximize efficiency, we choose the largest possible block size L.
L = N - L_h + 1
print(f"1. The block size (L) for the long sequence is calculated as:")
print(f"   L = N - L_h + 1 = {N} - {L_h} + 1 = {L}")

# The number of blocks is the total length of the long sequence divided by the block size.
# We use ceiling to ensure the entire sequence is covered.
num_blocks_add = math.ceil(L_x / L)
print(f"2. The number of non-overlapping blocks (K_add) needed is:")
print(f"   K_add = ceil(L_x / L) = ceil({L_x} / {L}) = {num_blocks_add}")

print(f"\nFor the overlap-add method, {num_blocks_add} blocks are processed.")
print(f"This requires {num_blocks_add} (DFT+IDFT) operations.")
print("-" * 30)

# Step 3: Calculate operations for Overlap-Save method
print("--- Overlap-Save Method Calculation ---")
# In this method, each processed block of length N produces M valid output samples.
# The number of valid samples M is determined by the same relationship.
M = N - L_h + 1
print(f"1. The number of valid output samples per block (M) is:")
print(f"   M = N - L_h + 1 = {N} - {L_h} + 1 = {M}")

# We need to generate the complete linear convolution output.
# The total length of the linear convolution is L_y = L_x + L_h - 1.
L_y = L_x + L_h - 1
print(f"2. The total length of the convolution output (L_y) is:")
print(f"   L_y = L_x + L_h - 1 = {L_x} + {L_h} - 1 = {L_y}")

# The number of blocks is the total output length divided by the valid samples per block.
# We use ceiling to ensure the entire output sequence is generated.
num_blocks_save = math.ceil(L_y / M)
print(f"3. The number of overlapping blocks (K_save) needed is:")
print(f"   K_save = ceil(L_y / M) = ceil({L_y} / {M}) = {num_blocks_save}")

print(f"\nFor the overlap-save method, {num_blocks_save} blocks are processed to get the complete output.")
print(f"This requires {num_blocks_save} (DFT+IDFT) operations.")
print("-" * 30)

# Step 4: Final Summary
print("\n--- Summary ---")
print(f"Overlap-Add requires: {num_blocks_add} (DFT+IDFT) operations.")
print(f"Overlap-Save requires: {num_blocks_save} (DFT+IDFT) operations.")

import math

# Step 1: Define the given parameters from the problem description.
# L is the length of the long sequence (x[n]).
# M is the length of the shorter sequence or filter (h[n]).
# N is the size of the DFT and IDFT used for the convolution.
L = 1200
M = 90
N = 128

print("--- Calculating required (DFT+IDFT) operations ---")

# Step 2: Calculate the number of operations for the Overlap-Add method.
# For overlap-add, the input signal is divided into non-overlapping blocks.
# The size of each data block (L_block) is chosen to be as large as possible
# such that the linear convolution result fits within the DFT size N.
# The condition is: L_block + M - 1 <= N, so we choose L_block = N - M + 1.
L_block = N - M + 1

# The number of operations is the total input length L divided by the data block length,
# rounded up to the nearest integer (ceiling), as even a partial last block needs a full operation.
num_ops_add = math.ceil(L / L_block)

print("\nOverlap-Add Method:")
print(f"The number of operations is calculated as ceil(L / (N - M + 1))")
print(f"Number of operations = ceil({L} / ({N} - {M} + 1)) = ceil({L} / {L_block}) = {num_ops_add}")

# Step 3: Calculate the number of operations for the Overlap-Save method.
# For overlap-save, we calculate the total length of the final output sequence, P.
# The output length of a linear convolution is P = L + M - 1.
P = L + M - 1

# Each (DFT+IDFT) operation produces a block of L_block valid (non-aliased) output samples.
# The size of this valid block is also N - M + 1.
# The number of operations is the total output length P divided by the valid samples per block,
# rounded up to the nearest integer (ceiling).
num_ops_save = math.ceil(P / L_block)

print("\nOverlap-Save Method:")
print(f"The number of operations is calculated as ceil((L + M - 1) / (N - M + 1))")
print(f"Total output length P = {L} + {M} - 1 = {P}")
print(f"Number of operations = ceil({P} / ({N} - {M} + 1)) = ceil({P} / {L_block}) = {num_ops_save}")
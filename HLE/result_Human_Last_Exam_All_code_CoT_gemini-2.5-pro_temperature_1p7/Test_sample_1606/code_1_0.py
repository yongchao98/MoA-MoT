import math

# Step 1: Define the parameters based on the problem description.
# Block convolution is used to filter a long signal with a shorter one.
# L_x is the length of the long signal, L_h is the length of the shorter filter.
L_x = 1200
L_h = 90

# N is the size of the DFT/IDFT used.
N = 128

# A necessary condition for block convolution is N >= L_h.
# 128 >= 90, so this condition is met.

# Step 2: Calculate the number of operations for the Overlap-Add method.
print("--- Overlap-Add Method Calculation ---")

# For Overlap-Add, we segment the input signal x[n] into non-overlapping
# blocks of length L_s. To avoid aliasing, the DFT size N must be >= L_s + L_h - 1.
# To minimize operations, we choose the maximum possible block size L_s.
L_s = N - L_h + 1
print(f"The optimal input block size L_s = N - L_h + 1 = {N} - {L_h} + 1 = {L_s}")

# The number of blocks (K_add) determines the number of (DFT+IDFT) operations.
# It is calculated by dividing the input signal length by the block size and rounding up.
K_add = math.ceil(L_x / L_s)
print("The number of required (DFT+IDFT) operations is K_add = ceil(L_x / L_s)")
print(f"K_add = ceil({L_x} / {L_s}) = {K_add}")
print("-" * 38 + "\n")


# Step 3: Calculate the number of operations for the Overlap-Save method.
print("--- Overlap-Save Method Calculation ---")

# For Overlap-Save, each N-point block processing yields L_s = N - L_h + 1 valid output samples.
# The value of L_s is the same as calculated before.
print(f"The number of valid output samples per block is L_s = N - L_h + 1 = {N} - {L_h} + 1 = {L_s}")

# We need to generate a total output sequence of length L_y.
L_y = L_x + L_h - 1
print(f"The total length of the output signal is L_y = L_x + L_h - 1 = {L_x} + {L_h} - 1 = {L_y}")

# The number of blocks (K_save) is the number required to generate the full output.
# It is calculated by dividing the total output length by the valid samples per block and rounding up.
K_save = math.ceil(L_y / L_s)
print("The number of required (DFT+IDFT) operations is K_save = ceil(L_y / L_s)")
print(f"K_save = ceil({L_y} / {L_s}) = {K_save}")
print("-" * 38)
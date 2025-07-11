import math

#
# This script calculates the number of (DFT+IDFT) operations required for
# fast convolution using the Overlap-Add and Overlap-Save methods.
#

# --- Given Parameters ---
L = 1200  # Length of the long sequence x[n]
M = 90    # Length of the shorter sequence h[n] (the filter)
N = 128   # The size of the DFT and IDFT used for the computation

# --- Calculation for both Overlap-Add and Overlap-Save ---

# For both methods, the number of new or useful data points (L_data)
# that can be processed in each block is determined by the constraint
# that the DFT size N must be large enough to hold the result of the
# linear convolution of a data block with the filter without aliasing.
# The length of the linear convolution of a block of length L_data and a
# filter of length M is L_data + M - 1.
# So, the constraint is N >= L_data + M - 1.
# To maximize efficiency, we choose the largest possible L_data:
# L_data = N - M + 1

L_data = N - M + 1

# Check if the DFT size is valid. If L_data is not positive, the DFT size N
# is too small for the given filter length M.
if L_data <= 0:
    print(f"Error: The DFT size N={N} is too small for the filter of length M={M}.")
    print(f"The condition N >= M must be met for overlap-save, and N >= M-1+L_data must be met for overlap-add with L_data > 0.")
else:
    # The total number of blocks required is the total length of the signal (L)
    # divided by the number of new data points processed per block (L_data).
    # We use math.ceil to round up to ensure all data is processed.
    num_blocks = math.ceil(L / L_data)

    # Each block requires one DFT and one IDFT. The DFT of the filter h[n] is
    # pre-computed once and reused, so it's not part of the per-block operations.
    # Therefore, the number of (DFT+IDFT) operations is equal to the number of blocks.

    # --- Print the results ---
    print("For an overlap-add implementation:")
    # We explicitly show the numbers in the final equation.
    print(f"Number of (DFT+IDFT) operations = ceil(L / (N - M + 1)) = ceil({L} / ({N} - {M} + 1)) = ceil({L} / {L_data}) = {num_blocks}")

    print("\nFor an overlap-save implementation:")
    # The calculation is identical for overlap-save.
    print(f"Number of (DFT+IDFT) operations = ceil(L / (N - M + 1)) = ceil({L} / ({N} - {M} + 1)) = ceil({L} / {L_data}) = {num_blocks}")

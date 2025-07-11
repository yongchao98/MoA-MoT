import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods.
    """
    # Given parameters
    # Let the shorter sequence be the filter h[n] and the longer be the signal x[n]
    L = 90      # Length of the shorter sequence h[n]
    M = 1200    # Length of the longer sequence x[n]
    N = 128     # DFT/IDFT size

    # --- Overlap-Add Method ---
    # For overlap-add, the input signal x[n] is segmented into non-overlapping
    # blocks of length L_block. To avoid aliasing, the DFT size N must be
    # greater than or equal to the convolution result length: N >= L_block + L - 1.
    # To maximize efficiency, we use the largest possible block length.
    l_block_add = N - L + 1

    # The number of blocks determines the number of (DFT+IDFT) operations.
    # We need to cover the entire length M of the signal x[n].
    num_ops_add = math.ceil(M / l_block_add)

    print("Overlap-Add Method Calculation:")
    print(f"Length of sequence x[n], M = {M}")
    print(f"Length of sequence h[n], L = {L}")
    print(f"DFT size, N = {N}")
    print(f"Maximum block length, L_block = N - L + 1 = {N} - {L} + 1 = {l_block_add}")
    print(f"Number of (DFT+IDFT) operations = ceil(M / L_block) = ceil({M} / {l_block_add}) = {num_ops_add}")
    print("-" * 40)

    # --- Overlap-Save Method ---
    # For overlap-save, the input signal x[n] is segmented into overlapping
    # blocks of length N. Each circular convolution yields a number of valid
    # (non-aliased) samples.
    # Number of valid samples per block = N - L + 1.
    l_valid_save = N - L + 1

    # The number of blocks needed is determined by how many blocks are required
    # to produce at least M valid output samples.
    num_ops_save = math.ceil(M / l_valid_save)

    print("Overlap-Save Method Calculation:")
    print(f"Length of sequence x[n], M = {M}")
    print(f"Length of sequence h[n], L = {L}")
    print(f"DFT size, N = {N}")
    print(f"Valid points per block, L_valid = N - L + 1 = {N} - {L} + 1 = {l_valid_save}")
    print(f"Number of (DFT+IDFT) operations = ceil(M / L_valid) = ceil({M} / {l_valid_save}) = {num_ops_save}")

solve_convolution_operations()
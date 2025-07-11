import math

def calculate_fft_operations():
    """
    Calculates the number of DFT and IDFT operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Parameters from the problem statement
    len_seq1 = 90
    len_seq2 = 1200
    N = 128  # DFT/IDFT size

    # For overlap-add/save, the filter must be shorter than the DFT length.
    # Therefore, the shorter sequence is the filter (M) and the longer
    # is the input signal (L).
    M = min(len_seq1, len_seq2)
    L = max(len_seq1, len_seq2)

    print(f"Problem Parameters:")
    print(f"Input Signal Length (L): {L}")
    print(f"Filter Length (M): {M}")
    print(f"DFT/IDFT Size (N): {N}")
    print("-" * 30)

    # Both methods process a certain number of new samples per block.
    # This value, L_block, is determined by N and M.
    # L_block = N - M + 1
    l_block = N - M + 1

    # The number of blocks (K) is determined by how many blocks of size L_block
    # are needed to process the entire input signal of length L.
    K = math.ceil(L / l_block)

    # --- Overlap-Add Calculation ---
    # Total ops = 1 (for filter DFT) + K (for input blocks DFTs) + K (for output blocks IDFTs)
    total_ops_add = 1 + 2 * K

    print("Overlap-Add Method:")
    print(f"  Number of new samples per block (L_block = N - M + 1): {N} - {M} + 1 = {l_block}")
    print(f"  Number of blocks required (K = ceil(L / L_block)): ceil({L} / {l_block}) = {K}")
    print(f"  Total operations = 1 (filter DFT) + {K} (input DFTs) + {K} (output IDFTs)")
    print(f"  Total = 1 + {2 * K} = {total_ops_add}")
    print("-" * 30)


    # --- Overlap-Save Calculation ---
    # The logic is identical for calculating the number of operations.
    # Each block processes L_block = N - M + 1 new samples.
    # The number of blocks K is the same.
    # Total ops = 1 (for filter DFT) + K (for input blocks DFTs) + K (for output blocks IDFTs)
    total_ops_save = 1 + 2 * K

    print("Overlap-Save Method:")
    print(f"  Number of new samples per block (N - M + 1): {N} - {M} + 1 = {l_block}")
    print(f"  Number of blocks required (K = ceil(L / L_block)): ceil({L} / {l_block}) = {K}")
    print(f"  Total operations = 1 (filter DFT) + {K} (input DFTs) + {K} (output IDFTs)")
    print(f"  Total = 1 + {2 * K} = {total_ops_save}")

# Execute the function to print the results
calculate_fft_operations()
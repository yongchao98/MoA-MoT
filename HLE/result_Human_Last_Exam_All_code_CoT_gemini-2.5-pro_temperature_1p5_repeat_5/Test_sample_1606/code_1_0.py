import math

def calculate_dft_operations():
    """
    Calculates the number of (DFT+IDFT) operations required for linear convolution
    using overlap-add and overlap-save methods.
    """
    len_seq1 = 90
    len_seq2 = 1200
    N_dft = 128

    # The filter (h[n]) is the shorter sequence, its length is M.
    # The signal (x[n]) is the longer sequence, its length is L_x.
    M = min(len_seq1, len_seq2)
    L_x = max(len_seq1, len_seq2)

    # Check if DFT size N is valid. We must have N >= M.
    if N_dft < M:
        print(f"Error: DFT size N={N_dft} must be greater than or equal to the filter length M={M}.")
        return

    print(f"Given Parameters:")
    print(f"Signal Length (L_x): {L_x}")
    print(f"Filter Length (M): {M}")
    print(f"DFT Size (N): {N_dft}")
    print("-" * 30)

    # --- Overlap-Add Method ---
    # The number of non-overlapping samples in each input block.
    L = N_dft - M + 1
    
    # The number of blocks needed is ceil(L_x / L).
    # Each block requires one DFT and one IDFT operation.
    num_ops_add = math.ceil(L_x / L)

    print("Overlap-Add Method:")
    print(f"The length of each non-overlapping data block (L) is calculated as:")
    print(f"L = N - M + 1 = {N_dft} - {M} + 1 = {L}")
    print(f"\nThe number of required (DFT+IDFT) operations is:")
    print(f"Ops_add = ceil(L_x / L) = ceil({L_x} / {L}) = {num_ops_add}")
    print("-" * 30)

    # --- Overlap-Save Method ---
    # The number of valid output samples from each block.
    # The formula for L is the same as in overlap-add.
    L = N_dft - M + 1
    
    # The total length of the output convolution result.
    L_y = L_x + M - 1
    
    # The number of blocks needed is ceil(L_y / L).
    # Each block requires one DFT and one IDFT operation.
    num_ops_save = math.ceil(L_y / L)
    
    print("Overlap-Save Method:")
    print(f"The number of valid output samples from each block (L) is:")
    print(f"L = N - M + 1 = {N_dft} - {M} + 1 = {L}")
    print(f"\nThe total length of the output sequence (L_y) is:")
    print(f"L_y = L_x + M - 1 = {L_x} + {M} - 1 = {L_y}")
    print(f"\nThe number of required (DFT+IDFT) operations is:")
    print(f"Ops_save = ceil(L_y / L) = ceil({L_y} / {L}) = {num_ops_save}")

calculate_dft_operations()
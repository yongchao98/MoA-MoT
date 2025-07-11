import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for convolving two sequences.
    """
    # 1. Identify parameters
    # The longer sequence is treated as the input signal 'x' to be segmented.
    # The shorter sequence is treated as the filter impulse response 'h'.
    L_x = 1200  # Length of the long sequence (input signal)
    M = 90      # Length of the short sequence (filter)
    N = 128     # DFT/IDFT size

    print("--- Problem Parameters ---")
    print(f"Length of the long sequence to be processed (L_x): {L_x}")
    print(f"Length of the short sequence (M): {M}")
    print(f"DFT size (N): {N}")
    print("-" * 30)

    # For both methods, the number of new/non-overlapping data points processed per block is L.
    # This value is constrained by the DFT size N and the filter length M.
    # L = N - M + 1
    L = N - M + 1
    
    # --- Overlap-Add Method Calculation ---
    print("--- Overlap-Add Method ---")
    # In overlap-add, the input signal is segmented into K non-overlapping blocks of length L.
    # Each block is then padded to length N for the DFT.
    # The number of blocks K determines the number of (DFT+IDFT) operations.
    # K = ceil(Total length / Length of each block)
    num_blocks_add = math.ceil(L_x / L)

    print(f"The size of each non-overlapping data block (L) is N - M + 1.")
    print(f"Equation: L = {N} - {M} + 1 = {L}")
    print(f"The number of blocks needed is ceil(L_x / L).")
    print(f"Equation: Number of operations = ceil({L_x} / {L}) = {num_blocks_add}")
    print("-" * 30)

    # --- Overlap-Save Method Calculation ---
    print("--- Overlap-Save Method ---")
    # In overlap-save, the input is segmented into K overlapping blocks of length N.
    # Each block produces L valid output samples, corresponding to L new input samples.
    # The number of blocks K required to process the whole signal is K = ceil(L_x / L).
    num_blocks_save = math.ceil(L_x / L)

    print(f"The number of new data points per block (L) is N - M + 1.")
    print(f"Equation: L = {N} - {M} + 1 = {L}")
    print(f"The number of blocks needed is ceil(L_x / L).")
    print(f"Equation: Number of operations = ceil({L_x} / {L}) = {num_blocks_save}")
    print("-" * 30)

    print(f"Final Answer:")
    print(f"Overlap-add operations: {num_blocks_add}")
    print(f"Overlap-save operations: {num_blocks_save}")

solve_convolution_operations()
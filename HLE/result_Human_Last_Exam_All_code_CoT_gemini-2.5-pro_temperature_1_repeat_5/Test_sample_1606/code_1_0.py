import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Given parameters
    L = 1200  # Length of the long sequence x[n]
    M = 90    # Length of the short sequence h[n]
    N = 128   # DFT/IDFT size

    print("Given Parameters:")
    print(f"Length of long sequence (L): {L}")
    print(f"Length of short sequence (M): {M}")
    print(f"DFT/IDFT size (N): {N}\n")

    # --- Overlap-Add Method Calculation ---
    print("--- Overlap-Add Implementation ---")
    
    # In overlap-add, the input sequence is broken into non-overlapping blocks of length L_b.
    # To avoid aliasing, the DFT size N must be at least L_b + M - 1.
    # We choose the maximum possible block size L_b = N - M + 1.
    L_b = N - M + 1
    print(f"The length of each data block is calculated as L_b = N - M + 1.")
    print(f"L_b = {N} - {M} + 1 = {L_b}")

    # The number of blocks K_add is the total length L divided by L_b, rounded up.
    K_add_float = L / L_b
    K_add = math.ceil(K_add_float)
    print(f"\nThe number of blocks K_add is calculated by ceil(L / L_b).")
    print(f"K_add = ceil({L} / {L_b}) = ceil({K_add_float:.3f}) = {K_add}")
    print(f"Total (DFT+IDFT) operations for Overlap-Add: {K_add}\n")


    # --- Overlap-Save Method Calculation ---
    print("--- Overlap-Save Implementation ---")
    
    # In overlap-save, the input is broken into overlapping blocks of length N.
    # Each block produces N - M + 1 valid output samples. This is also the
    # number of new samples from the input sequence processed in each step after the first.
    L_valid = N - M + 1
    print(f"The number of valid (new) samples processed per block is N - M + 1.")
    print(f"L_valid = {N} - {M} + 1 = {L_valid}")

    # We need to find the number of blocks K_save to process the entire sequence of length L.
    # The total number of samples processed by K blocks is N (for the first block)
    # plus (K - 1) * L_valid for the subsequent blocks. This must be >= L.
    # N + (K_save - 1) * L_valid >= L
    L_minus_N = L - N
    K_save_inequality_rhs = L_minus_N / L_valid
    K_save_final_rhs = K_save_inequality_rhs + 1
    K_save = math.ceil(K_save_final_rhs)

    print("\nTo find the number of blocks K_save, we solve the inequality:")
    print(f"N + (K_save - 1) * L_valid >= L")
    print(f"{N} + (K_save - 1) * {L_valid} >= {L}")
    print(f"(K_save - 1) * {L_valid} >= {L} - {N}")
    print(f"K_save - 1 >= {L_minus_N} / {L_valid}")
    print(f"K_save - 1 >= {K_save_inequality_rhs:.3f}")
    print(f"K_save >= {K_save_final_rhs:.3f}")
    print(f"The smallest integer K_save is {K_save}.")
    print(f"Total (DFT+IDFT) operations for Overlap-Save: {K_save}")

solve_convolution_operations()
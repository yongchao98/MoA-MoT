import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Given parameters
    L = 90    # Length of the first sequence
    M = 1200  # Length of the second (long) sequence
    N = 128   # DFT size

    print("Problem Parameters:")
    print(f"Length of shorter sequence, L = {L}")
    print(f"Length of longer sequence, M = {M}")
    print(f"DFT/IDFT size, N = {N}\n")

    # --- Overlap-Add Method Calculation ---
    print("--- Overlap-Add Method ---")
    
    # To avoid aliasing, block size L_b must satisfy: L_b + L - 1 <= N
    # To minimize blocks, we use the maximum possible block size.
    L_b = N - L + 1
    print(f"Maximum block size (L_b) = N - L + 1 = {N} - {L} + 1 = {L_b}")

    if L_b <= 0:
        print("Error: DFT size N is too small for overlap-add.")
        K_add = float('inf')
    else:
        # Number of blocks is ceil(M / L_b)
        K_add = math.ceil(M / L_b)
        print(f"Number of blocks (K_add) = ceil(M / L_b) = ceil({M} / {L_b}) = {K_add}")
    
    print(f"Total (DFT+IDFT) operations for Overlap-Add: {K_add}\n")

    # --- Overlap-Save Method Calculation ---
    print("--- Overlap-Save Method ---")
    
    # Each block produces L_useful valid output samples
    L_useful = N - L + 1
    print(f"Useful samples per block (L_useful) = N - L + 1 = {N} - {L} + 1 = {L_useful}")
    
    # Total length of the output sequence
    L_out = M + L - 1
    print(f"Total output length required (L_out) = M + L - 1 = {M} + {L} - 1 = {L_out}")

    if L_useful <= 0:
        print("Error: DFT size N is too small for overlap-save.")
        K_save = float('inf')
    else:
        # Number of blocks is ceil(L_out / L_useful)
        K_save = math.ceil(L_out / L_useful)
        print(f"Number of blocks (K_save) = ceil(L_out / L_useful) = ceil({L_out} / {L_useful}) = {K_save}")
        
    print(f"Total (DFT+IDFT) operations for Overlap-Save: {K_save}\n")

solve_convolution_operations()
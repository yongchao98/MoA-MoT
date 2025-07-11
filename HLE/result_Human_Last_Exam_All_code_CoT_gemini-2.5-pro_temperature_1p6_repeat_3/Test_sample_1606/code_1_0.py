import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Parameters
    M = 90       # Length of the shorter sequence (filter)
    Nx = 1200    # Length of the longer sequence (signal)
    N_dft = 128  # Size of the DFT/IDFT

    print(f"Problem Parameters:")
    print(f"  Length of filter sequence (M) = {M}")
    print(f"  Length of signal sequence (Nx) = {Nx}")
    print(f"  DFT/IDFT size (N) = {N_dft}\n")

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method ---")
    # For overlap-add, the input signal is segmented into non-overlapping blocks
    # of length L. To avoid aliasing, L is chosen such that N >= L + M - 1.
    # We choose the maximum possible L.
    L_add = N_dft - M + 1
    print(f"The length of each non-overlapping data block (L) is calculated as:")
    print(f"  L = N - M + 1 = {N_dft} - {M} + 1 = {L_add}")

    # The number of blocks (K_add) determines the number of DFT/IDFT operations.
    # We use math.ceil because even a partial final block requires a full operation.
    K_add = math.ceil(Nx / L_add)
    print(f"The number of (DFT+IDFT) operations is the number of blocks needed to process the signal:")
    print(f"  Operations = ceil(Nx / L) = ceil({Nx} / {L_add}) = {K_add}\n")

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method ---")
    # For overlap-save, the input signal is segmented into overlapping blocks
    # of length N. Each block produces L valid output samples.
    L_save = N_dft - M + 1
    print(f"The number of valid output samples from each block (L) is calculated as:")
    print(f"  L = N - M + 1 = {N_dft} - {M} + 1 = {L_save}")

    # The number of blocks (K_save) is determined by how many blocks are needed.
    # Each block effectively processes L new samples.
    K_save = math.ceil(Nx / L_save)
    print(f"The number of (DFT+IDFT) operations is the number of blocks needed to process the signal:")
    print(f"  Operations = ceil(Nx / L) = ceil({Nx} / {L_save}) = {K_save}")

solve_convolution_operations()
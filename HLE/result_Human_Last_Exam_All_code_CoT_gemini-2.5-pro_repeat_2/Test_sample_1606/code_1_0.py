import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Given parameters from the problem
    Lx = 1200  # Length of the long sequence x[n]
    M = 90     # Length of the short sequence h[n]
    N = 128    # DFT/IDFT size

    print("--- Calculating Operations for Overlap-Add Method ---")

    # For overlap-add, the input signal is broken into non-overlapping blocks of length L.
    # The condition N >= L + M - 1 must hold for the circular convolution
    # to be equivalent to the linear convolution of the block.
    # We choose the maximum possible block length L to maximize efficiency.
    L_add = N - M + 1
    print(f"1. The size of each non-overlapping input block (L) is calculated as N - M + 1.")
    print(f"   L = {N} - {M} + 1 = {L_add}")

    # The number of blocks is the total input length divided by the block length, rounded up.
    K_add = math.ceil(Lx / L_add)
    print(f"2. The number of required operations (K_add) is ceil(Lx / L).")
    print(f"   K_add = ceil({Lx} / {L_add}) = {K_add}")
    print(f"\nResult for Overlap-Add: {K_add} (DFT+IDFT) operations are necessary.\n")


    print("--- Calculating Operations for Overlap-Save Method ---")

    # For overlap-save, each block of length N produces L = N - M + 1 valid output samples.
    L_save = N - M + 1
    print(f"1. The number of valid output samples per block (L) is calculated as N - M + 1.")
    print(f"   L = {N} - {M} + 1 = {L_save}")


    # First, we determine the total length of the output sequence from linear convolution.
    Ly = Lx + M - 1
    print(f"2. The total length of the final output sequence (Ly) is Lx + M - 1.")
    print(f"   Ly = {Lx} + {M} - 1 = {Ly}")

    # The number of blocks is the total output length divided by the valid samples per block, rounded up.
    K_save = math.ceil(Ly / L_save)
    print(f"3. The number of required operations (K_save) is ceil(Ly / L).")
    print(f"   K_save = ceil({Ly} / {L_save}) = {K_save}")
    print(f"\nResult for Overlap-Save: {K_save} (DFT+IDFT) operations are necessary.")

if __name__ == "__main__":
    solve_convolution_operations()
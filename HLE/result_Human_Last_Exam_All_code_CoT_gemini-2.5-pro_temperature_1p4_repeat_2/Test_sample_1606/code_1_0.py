import math

def calculate_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # --- Given Parameters ---
    N1 = 90
    N2 = 1200
    N_dft = 128

    # Identify the long sequence (x[n]) and the short sequence (h[n], the filter).
    # Convolution is commutative, so we break up the longer sequence for efficiency.
    L_x = max(N1, N2)
    M = min(N1, N2)

    # --- Overlap-Add Method ---
    # Determine the data block length L that avoids time-domain aliasing.
    # N_dft >= L + M - 1  =>  L <= N_dft - M + 1
    L_block_add = N_dft - M + 1

    # The number of blocks (K) is the total length of the long sequence
    # divided by the block length, rounded up (ceiling).
    K_add = math.ceil(L_x / L_block_add)

    # --- Overlap-Save Method ---
    # Determine the number of useful (non-aliased) samples per processed block.
    L_useful_save = N_dft - M + 1

    # The number of blocks (K) is the total length of the long sequence
    # divided by the number of useful samples per block, rounded up.
    K_save = math.ceil(L_x / L_useful_save)


    # --- Print the results step-by-step ---
    print("--- Overlap-Add Implementation ---")
    print(f"Length of long sequence (L_x) = {L_x}")
    print(f"Length of short sequence (M) = {M}")
    print(f"DFT size (N) = {N_dft}")
    print("\nCalculation:")
    print(f"1. Data block length (L) = N - M + 1 = {N_dft} - {M} + 1 = {L_block_add}")
    print(f"2. Number of (DFT+IDFT) operations = ceil(L_x / L) = ceil({L_x} / {L_block_add}) = {K_add}")
    print(f"\nResult: Overlap-add requires {K_add} (DFT+IDFT) operations.")
    print("-" * 40)

    print("\n--- Overlap-Save Implementation ---")
    print(f"Length of long sequence (L_x) = {L_x}")
    print(f"Length of short sequence (M) = {M}")
    print(f"DFT size (N) = {N_dft}")
    print("\nCalculation:")
    print(f"1. Useful samples per block (L_useful) = N - M + 1 = {N_dft} - {M} + 1 = {L_useful_save}")
    print(f"2. Number of (DFT+IDFT) operations = ceil(L_x / L_useful) = ceil({L_x} / {L_useful_save}) = {K_save}")
    print(f"\nResult: Overlap-save requires {K_save} (DFT+IDFT) operations.")
    print("-" * 40)


if __name__ == '__main__':
    calculate_convolution_operations()

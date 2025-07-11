import math

def calculate_fft_operations():
    """
    Calculates the number of (DFT+IDFT) operations for Overlap-Add and
    Overlap-Save methods for linear convolution.
    """
    # Parameters from the problem statement
    # Let L_x be the length of the long sequence and L_h be the length of the short sequence.
    L_x = 1200
    L_h = 90
    N = 128 # N-point DFT/IDFT

    print(f"Given parameters:")
    print(f"  Length of long sequence (L_x): {L_x}")
    print(f"  Length of short sequence (L_h): {L_h}")
    print(f"  DFT/IDFT size (N): {N}\n")

    # --- Overlap-Add Method ---
    print("Overlap-Add Method Calculation:")
    # For Overlap-Add, we divide the long sequence into non-overlapping blocks of length L_b.
    # The linear convolution of a block of length L_b and a filter of length L_h
    # results in an output of length L_b + L_h - 1. This must be <= N.
    # To maximize efficiency, we choose the largest possible L_b.
    # L_b + L_h - 1 = N => L_b = N - L_h + 1
    L_b_add = N - L_h + 1
    print(f"  Length of each data block (L_b) = N - L_h + 1 = {N} - {L_h} + 1 = {L_b_add}")

    # Number of blocks needed is the ceiling of the total length divided by block length.
    K_add = math.ceil(L_x / L_b_add)
    print(f"  Number of blocks (K) = ceil(L_x / L_b) = ceil({L_x} / {L_b_add}) = {K_add}")

    # Total operations = 1 DFT for the filter + K DFTs for data blocks + K IDFTs for results.
    dfts_add = 1 + K_add
    idfts_add = K_add
    total_ops_add = dfts_add + idfts_add
    print(f"  Total operations = (1 DFT for filter + {K_add} DFTs for data) + ({K_add} IDFTs for results)")
    print(f"  Total operations = {dfts_add} DFTs + {idfts_add} IDFTs = {total_ops_add}\n")

    # --- Overlap-Save Method ---
    print("Overlap-Save Method Calculation:")
    # For Overlap-Save, we use overlapping input blocks of length N.
    # Each block produces N output samples, but only L_b = N - L_h + 1 are valid (not aliased).
    L_b_save = N - L_h + 1
    print(f"  Number of new/valid samples per block (L_b) = N - L_h + 1 = {N} - {L_h} + 1 = {L_b_save}")

    # Number of blocks needed is based on how many sets of L_b new samples are required.
    K_save = math.ceil(L_x / L_b_save)
    print(f"  Number of blocks (K) = ceil(L_x / L_b) = ceil({L_x} / {L_b_save}) = {K_save}")

    # Total operations = 1 DFT for the filter + K DFTs for data blocks + K IDFTs for results.
    dfts_save = 1 + K_save
    idfts_save = K_save
    total_ops_save = dfts_save + idfts_save
    print(f"  Total operations = (1 DFT for filter + {K_save} DFTs for data) + ({K_save} IDFTs for results)")
    print(f"  Total operations = {dfts_save} DFTs + {idfts_save} IDFTs = {total_ops_save}")

calculate_fft_operations()
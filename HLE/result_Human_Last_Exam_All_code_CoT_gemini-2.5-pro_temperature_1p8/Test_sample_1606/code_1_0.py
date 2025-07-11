import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for Overlap-Add and Overlap-Save methods.
    """
    # Step 1: Define problem parameters
    # The standard convention is to filter a long sequence (L) with a shorter one (M).
    len_x = 1200  # Length of the long sequence, L
    len_h = 90    # Length of the short sequence, M
    dft_size = 128 # DFT/IDFT size, N

    # --- Overlap-Add Method Calculation ---
    
    # Step 2.1: Determine the processing block size (Lx) for non-overlapping blocks.
    # To prevent time-domain aliasing in the N-point circular convolution,
    # the block size Lx must be N - M + 1.
    block_size_add = dft_size - len_h + 1

    # Step 2.2: Calculate the number of blocks required. This corresponds to the number
    # of (DFT+IDFT) operations. We use ceiling to account for the last partial block.
    num_ops_add = math.ceil(len_x / block_size_add)


    # --- Overlap-Save Method Calculation ---

    # Step 3.1: Determine the number of useful/valid output samples per block (Lu).
    # For each N-length input block, N-M+1 samples are valid.
    useful_samples_save = dft_size - len_h + 1
    
    # Step 3.2: Calculate the number of blocks required. This corresponds to the number
    # of (DFT+IDFT) operations. Each block processes 'useful_samples_save' new samples from the input.
    num_ops_save = math.ceil(len_x / useful_samples_save)

    # --- Print Results ---
    print("--- Problem Parameters ---")
    print(f"Length of long sequence (L): {len_x}")
    print(f"Length of short sequence (M): {len_h}")
    print(f"DFT size (N): {dft_size}\n")

    print("--- Overlap-Add Method ---")
    print("The number of operations is ceil(L / (N - M + 1))")
    print(f"Calculation: ceil({len_x} / ({dft_size} - {len_h} + 1))")
    print(f"             = ceil({len_x} / {block_size_add})")
    print(f"             = {num_ops_add}")
    print(f"Result: {num_ops_add} (DFT+IDFT) operations are necessary for an overlap-add implementation.\n")

    print("--- Overlap-Save Method ---")
    print("The number of operations is ceil(L / (N - M + 1))")
    print(f"Calculation: ceil({len_x} / ({dft_size} - {len_h} + 1))")
    print(f"             = ceil({len_x} / {useful_samples_save})")
    print(f"             = {num_ops_save}")
    print(f"Result: {num_ops_save} (DFT+IDFT) operations are necessary for an overlap-save implementation.")

solve_convolution_operations()
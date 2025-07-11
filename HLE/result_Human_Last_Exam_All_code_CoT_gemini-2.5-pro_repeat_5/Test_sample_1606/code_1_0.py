import math

def calculate_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and overlap-save methods.
    """
    # Given parameters
    L_x = 1200  # Length of the long sequence x[n]
    L_h = 90    # Length of the short sequence h[n]
    N = 128     # DFT size

    print(f"Given parameters:")
    print(f"Length of long sequence (Lx): {L_x}")
    print(f"Length of short sequence (Lh): {L_h}")
    print(f"DFT size (N): {N}\n")

    # The DFT size N must be >= Lh for these methods to work.
    if N < L_h:
        print("Error: DFT size N must be greater than or equal to the filter length Lh.")
        return

    # Step 1: Calculate the number of new/valid samples per block
    L_block = N - L_h + 1
    print(f"Step 1: Calculate the number of valid samples per block (L_block)")
    print(f"L_block = N - Lh + 1 = {N} - {L_h} + 1 = {L_block}\n")

    # Step 2: Calculate operations for Overlap-Add method
    # The number of blocks is determined by sectioning the input sequence of length Lx.
    # Each block requires one (DFT+IDFT) operation.
    K_add = math.ceil(L_x / L_block)
    print("--- Overlap-Add Method ---")
    print("The number of operations is determined by the number of blocks needed to cover the input sequence.")
    print(f"Number of operations = ceil(Lx / L_block)")
    print(f"                   = ceil({L_x} / {L_block})")
    print(f"                   = ceil({L_x / L_block:.3f})")
    print(f"                   = {K_add}\n")

    # Step 3: Calculate operations for Overlap-Save method
    # The number of blocks is determined by sectioning the output sequence.
    # Total length of the output sequence Ly = Lx + Lh - 1
    L_y = L_x + L_h - 1
    # Each block requires one (DFT+IDFT) operation.
    K_save = math.ceil(L_y / L_block)
    print("--- Overlap-Save Method ---")
    print("The number of operations is determined by the number of blocks needed to generate the full output sequence.")
    print(f"Length of output sequence (Ly) = Lx + Lh - 1 = {L_x} + {L_h} - 1 = {L_y}")
    print(f"Number of operations = ceil(Ly / L_block)")
    print(f"                   = ceil({L_y} / {L_block})")
    print(f"                   = ceil({L_y / L_block:.3f})")
    print(f"                   = {K_save}\n")

    print("Final Answer:")
    print(f"Overlap-Add requires {K_add} (DFT+IDFT) operations.")
    print(f"Overlap-Save requires {K_save} (DFT+IDFT) operations.")

# Execute the function to print the results
calculate_convolution_operations()
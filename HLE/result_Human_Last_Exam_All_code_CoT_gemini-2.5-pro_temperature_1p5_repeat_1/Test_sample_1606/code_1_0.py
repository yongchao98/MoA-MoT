import math

def calculate_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for convolving two sequences.
    """
    # --- Given Parameters ---
    seq_len_1 = 90
    seq_len_2 = 1200
    N_dft = 128

    # Identify the long and short sequences
    if seq_len_1 > seq_len_2:
        L_long = seq_len_1
        L_short = seq_len_2
    else:
        L_long = seq_len_2
        L_short = seq_len_1

    print("--- Problem Setup ---")
    print(f"Length of the long sequence (L_long): {L_long}")
    print(f"Length of the short sequence (L_short): {L_short}")
    print(f"DFT size (N): {N_dft}\n")

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method Calculation ---")
    # Block size for the long sequence
    L_block_add = N_dft - L_short + 1
    print(f"1. Calculate block size = N - L_short + 1")
    print(f"   Block size = {N_dft} - {L_short} + 1 = {L_block_add}")
    
    # Number of operations for overlap-add
    num_ops_add = math.ceil(L_long / L_block_add)
    print(f"2. Calculate number of operations = ceil(L_long / Block size)")
    print(f"   Number of operations = ceil({L_long} / {L_block_add}) = {num_ops_add}")
    print(f"Result: Overlap-Add requires {num_ops_add} (DFT+IDFT) operations.\n")

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method Calculation ---")
    # Number of valid output samples per block
    L_valid_save = N_dft - L_short + 1
    print(f"1. Calculate valid samples per block = N - L_short + 1")
    print(f"   Valid samples = {N_dft} - {L_short} + 1 = {L_valid_save}")
    
    # Total length of the final convolved sequence
    L_output = L_long + L_short - 1
    print(f"2. Calculate total output length = L_long + L_short - 1")
    print(f"   Total output length = {L_long} + {L_short} - 1 = {L_output}")

    # Number of operations for overlap-save
    num_ops_save = math.ceil(L_output / L_valid_save)
    print(f"3. Calculate number of operations = ceil(Total output length / Valid samples)")
    print(f"   Number of operations = ceil({L_output} / {L_valid_save}) = {num_ops_save}")
    print(f"Result: Overlap-Save requires {num_ops_save} (DFT+IDFT) operations.\n")

if __name__ == "__main__":
    calculate_convolution_operations()

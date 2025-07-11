import math

def calculate_dft_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and overlap-save methods.
    """
    # Given parameters
    L = 90  # Length of the first sequence (e.g., filter h[n])
    M = 1200 # Length of the second sequence (e.g., signal x[n])
    N = 128  # DFT/IDFT size

    print(f"Given parameters:")
    print(f"Length of sequence h[n], L = {L}")
    print(f"Length of sequence x[n], M = {M}")
    print(f"DFT/IDFT size, N = {N}")
    print("-" * 30)

    # --- Overlap-Add Method ---
    print("Calculating for Overlap-Add Method:")
    # Block size P must satisfy: P + L - 1 <= N
    # So, the maximum block size is P = N - L + 1
    P_add = N - L + 1
    print(f"1. Maximum input block size (P) = N - L + 1 = {N} - {L} + 1 = {P_add}")

    # Number of blocks is ceil(M / P)
    num_blocks_add = math.ceil(M / P_add)
    print(f"2. Number of blocks = ceil(M / P) = ceil({M} / {P_add}) = {num_blocks_add}")
    print(f"For the overlap-add method, {num_blocks_add} (DFT+IDFT) operations are necessary.\n")

    # --- Overlap-Save Method ---
    print("Calculating for Overlap-Save Method:")
    # Number of valid output points per block
    P_valid_save = N - L + 1
    print(f"1. Valid points per block = N - L + 1 = {N} - {L} + 1 = {P_valid_save}")

    # Total length of the linear convolution output
    output_len = M + L - 1
    print(f"2. Total output length required = M + L - 1 = {M} + {L} - 1 = {output_len}")
    
    # Number of blocks is ceil(output_len / P_valid_save)
    num_blocks_save = math.ceil(output_len / P_valid_save)
    print(f"3. Number of blocks = ceil(Total output length / Valid points) = ceil({output_len} / {P_valid_save}) = {num_blocks_save}")
    print(f"For the overlap-save method, {num_blocks_save} (DFT+IDFT) operations are necessary.")
    print("-" * 30)

    # Final Answer Summary
    print(f"Final Answer:")
    print(f"Overlap-Add requires: {num_blocks_add} operations")
    print(f"Overlap-Save requires: {num_blocks_save} operations")

    # This is the final answer string requested by the prompt format
    global final_answer_string
    final_answer_string = f"<<<Overlap-Add: {num_blocks_add}, Overlap-Save: {num_blocks_save}>>>"


# Execute the calculation
calculate_dft_operations()
# The final answer is captured in a global variable and printed at the end
# to adhere to the prompt's output format.
print(final_answer_string)
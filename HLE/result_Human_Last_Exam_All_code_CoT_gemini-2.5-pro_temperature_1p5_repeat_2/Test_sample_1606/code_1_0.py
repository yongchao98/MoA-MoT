import math

def calculate_convolution_ops():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods.
    """
    # Given parameters
    N_h = 90   # Length of the shorter sequence (filter h[n]), M
    N_x = 1200 # Length of the longer sequence (signal x[n])
    N_dft = 128  # DFT size, N

    # For both methods, the number of useful/new samples per block is the same.
    # For overlap-add, this is the input block length.
    # For overlap-save, this is the valid output length.
    L = N_dft - N_h + 1

    # --- Overlap-Add Method ---
    # The number of operations is determined by how many blocks of length L
    # are needed to cover the entire sequence of length N_x.
    num_ops_add = math.ceil(N_x / L)

    print("Overlap-Add Method:")
    print(f"The input block length L is given by N_dft - N_h + 1.")
    print(f"L = {N_dft} - {N_h} + 1 = {L}")
    print(f"The number of (DFT+IDFT) operations is ceil(N_x / L).")
    print(f"Number of ops = ceil({N_x} / {L}) = {num_ops_add}")
    print("-" * 20)

    # --- Overlap-Save Method ---
    # The number of operations is determined by how many blocks are needed to
    # generate N_x useful samples, where each block generates L useful samples.
    num_ops_save = math.ceil(N_x / L)
    
    print("Overlap-Save Method:")
    print(f"The number of valid output samples L is given by N_dft - N_h + 1.")
    print(f"L = {N_dft} - {N_h} + 1 = {L}")
    print(f"The number of (DFT+IDFT) operations is ceil(N_x / L).")
    print(f"Number of ops = ceil({N_x} / {L}) = {num_ops_save}")
    print("-" * 20)
    
    # Although not explicitly asked for, we can form the final answer string here
    # global final_answer
    # final_answer = f"<<<{num_ops_add} for overlap-add and {num_ops_save} for overlap-save>>>"


# Run the calculation and print the results
calculate_convolution_ops()

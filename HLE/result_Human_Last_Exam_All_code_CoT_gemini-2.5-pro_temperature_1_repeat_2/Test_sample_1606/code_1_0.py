import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """

    # Given parameters
    L = 90      # Length of the shorter sequence (we'll treat this as the filter)
    M = 1200    # Length of the longer sequence (we'll treat this as the signal)
    N = 128     # DFT/IDFT size

    # --- Overlap-Add Method Calculation ---
    print("--- Overlap-Add Method ---")
    
    # In overlap-add, the input signal is divided into non-overlapping blocks of length K.
    # To prevent time-domain aliasing, K must satisfy: K + L - 1 <= N.
    # We choose the largest possible K.
    K_add = N - L + 1
    print(f"The maximum data block length (K) is determined by the formula: N - L + 1")
    print(f"K = {N} - {L} + 1 = {K_add}")

    # The number of blocks is the total signal length M divided by the block length K.
    # We use ceiling because a partial block still requires a full DFT/IDFT operation.
    num_blocks_add = math.ceil(M / K_add)
    print(f"The number of blocks is calculated by dividing the signal length M by K and taking the ceiling:")
    print(f"Number of blocks = ceil(M / K) = ceil({M} / {K_add}) = {num_blocks_add}")
    
    print(f"\nFor the Overlap-Add method, the number of (DFT+IDFT) operations is: {num_blocks_add}\n")

    # --- Overlap-Save Method Calculation ---
    print("--- Overlap-Save Method ---")

    # In overlap-save, each N-point circular convolution produces K valid output samples.
    K_save = N - L + 1
    print(f"The number of valid output samples per block (K) is given by: N - L + 1")
    print(f"K = {N} - {L} + 1 = {K_save}")

    # We need to generate the entire output sequence, whose length P is M + L - 1.
    P = M + L - 1
    print(f"The total length of the linear convolution output (P) is: M + L - 1")
    print(f"P = {M} + {L} - 1 = {P}")

    # The number of blocks needed is the total output length P divided by the valid samples per block K.
    # We use ceiling to ensure we generate all P samples.
    num_blocks_save = math.ceil(P / K_save)
    print(f"The number of blocks is calculated by dividing the total output length P by K and taking the ceiling:")
    print(f"Number of blocks = ceil(P / K) = ceil({P} / {K_save}) = {num_blocks_save}")
    
    print(f"\nFor the Overlap-Save method, the number of (DFT+IDFT) operations is: {num_blocks_save}")

solve_convolution_operations()
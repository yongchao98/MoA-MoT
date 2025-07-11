import math

def calculate_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for Overlap-Add and
    Overlap-Save methods for linear convolution.
    """
    # Given parameters
    L = 1200  # Length of the long sequence x[n]
    M = 90    # Length of the short sequence h[n]
    N = 128   # DFT/IDFT size

    print(f"Problem Parameters:")
    print(f"Length of sequence 'x[n]', L = {L}")
    print(f"Length of sequence 'h[n]', M = {M}")
    print(f"N-point DFT/IDFT size, N = {N}")
    print("-" * 30)

    # --- Overlap-Add Method Calculation ---
    # In the Overlap-Add method, the number of operations is determined by
    # partitioning the long input sequence L into non-overlapping blocks.
    print("1. Overlap-Add Method Calculation:")

    # Calculate the size of each input block
    l_block = N - M + 1
    print(f"   The input signal is split into blocks of size L_block.")
    print(f"   L_block = N - M + 1 = {N} - {M} + 1 = {l_block}")

    # Calculate the number of blocks needed to cover the input sequence
    k_add = math.ceil(L / l_block)
    print(f"   The number of operations is the number of blocks needed to cover the input signal.")
    print(f"   Number of operations = ceil(L / L_block) = ceil({L} / {l_block}) = {k_add}")
    print("-" * 30)

    # --- Overlap-Save Method Calculation ---
    # In the Overlap-Save method, the number of operations is determined by
    # the number of blocks required to generate the full output sequence.
    print("2. Overlap-Save Method Calculation:")

    # Calculate the length of the full linear convolution output
    l_y = L + M - 1
    print(f"   The total length of the linear convolution output is L_y.")
    print(f"   L_y = L + M - 1 = {L} + {M} - 1 = {l_y}")

    # Calculate the number of useful data points from each processed block
    l_data = N - M + 1
    print(f"   Each processed block yields L_data valid output samples.")
    print(f"   L_data = N - M + 1 = {N} - {M} + 1 = {l_data}")

    # Calculate the number of blocks needed to generate the full output
    k_save = math.ceil(l_y / l_data)
    print(f"   The number of operations is the number of blocks needed to generate the full output.")
    print(f"   Number of operations = ceil(L_y / L_data) = ceil({l_y} / {l_data}) = {k_save}")
    print("-" * 30)

    print("\nFinal Result:")
    print(f"Overlap-Add requires {k_add} (DFT+IDFT) operations.")
    print(f"Overlap-Save requires {k_save} (DFT+IDFT) operations.")


if __name__ == '__main__':
    calculate_convolution_operations()

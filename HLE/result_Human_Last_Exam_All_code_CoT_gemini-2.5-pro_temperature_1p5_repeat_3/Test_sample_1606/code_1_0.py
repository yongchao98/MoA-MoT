import math

def solve_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods.
    """
    # Given parameters
    L = 1200  # Length of the long sequence x[n]
    M = 90   # Length of the short sequence h[n]
    N = 128   # DFT/IDFT size

    # --- Overlap-Add Calculation ---
    # In overlap-add, we partition the input sequence x[n] of length L.
    # The number of new samples per block is N - M + 1.
    l_data_add = N - M + 1
    # Number of operations is the ceiling of L divided by the samples per block.
    k_add = math.ceil(L / l_data_add)

    print("Overlap-Add Method:")
    print(f"The number of useful data samples per block is N - M + 1 = {N} - {M} + 1 = {l_data_add}.")
    print("The number of operations is calculated by partitioning the input sequence of length L.")
    print(f"Number of operations = ceil(L / (N - M + 1)) = ceil({L} / {l_data_add}) = {k_add}")
    print("-" * 20)

    # --- Overlap-Save Calculation ---
    # In overlap-save, we construct the output sequence y[n].
    # The length of the output sequence y[n] is L + M - 1.
    l_y = L + M - 1
    # The number of valid output samples per block is N - M + 1.
    l_data_save = N - M + 1
    # Number of operations is the ceiling of the output length divided by valid samples per block.
    k_save = math.ceil(l_y / l_data_save)

    print("Overlap-Save Method:")
    print(f"The length of the output sequence is L + M - 1 = {L} + {M} - 1 = {l_y}.")
    print(f"The number of valid output samples per block is N - M + 1 = {N} - {M} + 1 = {l_data_save}.")
    print("The number of operations is calculated to generate the entire output sequence.")
    print(f"Number of operations = ceil((L + M - 1) / (N - M + 1)) = ceil({l_y} / {l_data_save}) = {k_save}")
    print("-" * 20)

    print(f"Final Answer: Overlap-Add requires {k_add} operations, and Overlap-Save requires {k_save} operations.")

solve_convolution_operations()
<<<Overlap-Add: 31, Overlap-Save: 34>>>
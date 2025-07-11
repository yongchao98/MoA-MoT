import math

def calculate_convolution_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Problem parameters
    len_seq1 = 90
    len_seq2 = 1200
    dft_size = 128

    # For overlap methods, the shorter sequence is the filter 'h' and the longer is the signal 'x'.
    M = min(len_seq1, len_seq2)
    L_x = max(len_seq1, len_seq2)
    N = dft_size

    # --- Overlap-Add Method ---
    # The signal x(n) is divided into non-overlapping blocks of length L_block.
    # L_block is chosen to maximize efficiency while satisfying N >= M + L_block - 1.
    # So, L_block = N - M + 1.
    # The number of operations is the total signal length divided by L_block, rounded up.
    
    L_block = N - M + 1
    num_ops_add = math.ceil(L_x / L_block)

    print("For an overlap-add implementation:")
    print(f"Number of (DFT+IDFT) operations = ceil(Signal Length / Block Length)")
    print(f"                               = ceil({L_x} / ({N} - {M} + 1))")
    print(f"                               = ceil({L_x} / {L_block})")
    print(f"                               = {num_ops_add}")

    # --- Overlap-Save Method ---
    # Each operation on an input block of length N yields L_valid = N - M + 1 useful output samples.
    # The total length of the final convolved sequence is L_y = L_x + M - 1.
    # The number of operations is the total output length divided by the number of valid
    # samples per operation, rounded up.

    L_y = L_x + M - 1
    L_valid = N - M + 1
    num_ops_save = math.ceil(L_y / L_valid)
    
    print("\nFor an overlap-save implementation:")
    print(f"Number of (DFT+IDFT) operations = ceil(Output Length / Valid Samples per Block)")
    print(f"                               = ceil(({L_x} + {M} - 1) / ({N} - {M} + 1))")
    print(f"                               = ceil({L_y} / {L_valid})")
    print(f"                               = {num_ops_save}")

if __name__ == '__main__':
    calculate_convolution_operations()

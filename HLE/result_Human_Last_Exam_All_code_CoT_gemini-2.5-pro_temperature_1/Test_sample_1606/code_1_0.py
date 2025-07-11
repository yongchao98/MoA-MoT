import math

def solve_convolution_operations():
    """
    Calculates the number of DFT and IDFT operations for overlap-add and
    overlap-save methods based on given sequence and DFT lengths.
    """
    # --- Given Parameters ---
    len_seq1 = 90
    len_seq2 = 1200
    N_dft = 128

    # By convention, the shorter sequence is the filter (M) and the longer
    # is the signal (Lx).
    M = min(len_seq1, len_seq2)
    Lx = max(len_seq1, len_seq2)

    print(f"Given Parameters:\n- Length of sequence 1: {len_seq1}\n- Length of sequence 2: {len_seq2}\n- DFT size (N): {N_dft}\n")
    print(f"Assuming shorter sequence is the filter (M={M}) and the longer is the signal (Lx={Lx}).\n")

    # --- Overlap-Add Method ---
    print("--- Overlap-Add Method Calculation ---")
    # Block length L_add is chosen such that L_add + M - 1 = N_dft
    L_add = N_dft - M + 1
    print(f"1. Data block length (L_add) = N - M + 1 = {N_dft} - {M} + 1 = {L_add}")

    # Number of blocks K_add is ceil(Lx / L_add)
    K_add = math.ceil(Lx / L_add)
    print(f"2. Number of blocks (K_add) = ceil(Lx / L_add) = ceil({Lx} / {L_add}) = {K_add}")

    # Total operations = 1 (for filter) + K_add (for data blocks) + K_add (for IDFTs)
    dfts_add = 1 + K_add
    idfts_add = K_add
    total_ops_add = dfts_add + idfts_add
    print(f"3. Total DFTs needed = 1 (for filter) + {K_add} (for blocks) = {dfts_add}")
    print(f"4. Total IDFTs needed = {K_add} (for blocks) = {idfts_add}")
    print(f"5. Total (DFT+IDFT) operations = {dfts_add} + {idfts_add} = {total_ops_add}\n")

    # --- Overlap-Save Method ---
    print("--- Overlap-Save Method Calculation ---")
    # Number of new data points per block L_save = N_dft - M + 1
    L_save = N_dft - M + 1
    print(f"1. New data points per block (L_save) = N - M + 1 = {N_dft} - {M} + 1 = {L_save}")

    # Number of blocks K_save is ceil(Lx / L_save)
    K_save = math.ceil(Lx / L_save)
    print(f"2. Number of blocks (K_save) = ceil(Lx / L_save) = ceil({Lx} / {L_save}) = {K_save}")

    # Total operations = 1 (for filter) + K_save (for data blocks) + K_save (for IDFTs)
    dfts_save = 1 + K_save
    idfts_save = K_save
    total_ops_save = dfts_save + idfts_save
    print(f"3. Total DFTs needed = 1 (for filter) + {K_save} (for blocks) = {dfts_save}")
    print(f"4. Total IDFTs needed = {K_save} (for blocks) = {idfts_save}")
    print(f"5. Total (DFT+IDFT) operations = {dfts_save} + {idfts_save} = {total_ops_save}\n")

    print(f"Final Answer: Overlap-Add requires {total_ops_add} operations, and Overlap-Save requires {total_ops_save} operations.")

solve_convolution_operations()
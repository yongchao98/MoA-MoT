import numpy as np

def calculate_transmitted_symbols():
    """
    This function calculates the transmitted symbols for a 2x2 MIMO system 
    with QPSK modulation and a modified Alamouti code.
    """
    # --- 1. Problem Setup ---
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j,
    }
    rotation_angle = np.pi / 8

    # --- 2. Bit to Symbol Mapping ---
    print("--- Step 1: Bit-to-Symbol Mapping ---")
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]
    
    print(f"Input bit stream: {bit_stream}")
    print(f"Grouped bit pairs: {bit_pairs}")
    print("Corresponding QPSK symbols (s_k):")
    for i, s in enumerate(qpsk_symbols):
        print(f"  s_{i+1} ({bit_pairs[i]}): {s}")

    # --- 3. Symbol Rotation ---
    print("\n--- Step 2: Symbol Rotation ---")
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]
    
    print(f"Each symbol is rotated by pi/8 radians (multiplied by {np.round(rotation_factor, 4)}).")
    print("Rotated symbols (x_k):")
    for i, x in enumerate(rotated_symbols):
        print(f"  x_{i+1} = {np.round(x, 4)}")
    
    # --- 4. Alamouti Encoding ---
    print("\n--- Step 3: Alamouti Encoding ---")
    print("Rotated symbols are processed in pairs (x_k, x_{k+1}) over two time slots.\n")
    tx_antenna_1 = []
    tx_antenna_2 = []

    for i in range(0, len(rotated_symbols), 2):
        x_k = rotated_symbols[i]
        x_k_plus_1 = rotated_symbols[i+1]
        
        # Populate symbol lists for antennas
        t1_t1 = x_k
        t2_t1 = x_k_plus_1
        t1_t2 = -np.conj(x_k_plus_1)
        t2_t2 = np.conj(x_k)
        
        tx_antenna_1.extend([t1_t1, t1_t2])
        tx_antenna_2.extend([t2_t1, t2_t2])
        
        # Print the equation for this pair
        print(f"Processing Pair {(i//2) + 1}: (x_{i+1}, x_{i+2})")
        print(f"  Time Slot {i+1}:")
        print(f"    Antenna 1 transmits x_{i+1} = {np.round(t1_t1, 4)}")
        print(f"    Antenna 2 transmits x_{i+2} = {np.round(t2_t1, 4)}")
        print(f"  Time Slot {i+2}:")
        print(f"    Antenna 1 transmits -x_{i+2}* = {np.round(t1_t2, 4)}")
        print(f"    Antenna 2 transmits x_{i+1}*  = {np.round(t2_t2, 4)}\n")
        
    # --- 5. Final Results ---
    print("\n--- Final Transmitted Symbol Sequences ---")
    print("Symbols from Antenna 1:")
    print("  " + ", ".join([f"({np.round(s, 4)})" for s in tx_antenna_1]))
    print("\nSymbols from Antenna 2:")
    print("  " + ", ".join([f"({np.round(s, 4)})" for s in tx_antenna_2]))


# Execute the calculation
calculate_transmitted_symbols()
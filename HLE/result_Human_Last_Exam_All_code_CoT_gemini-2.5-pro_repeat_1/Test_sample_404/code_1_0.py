import numpy as np

def solve_mimo_transmission():
    """
    Calculates and prints the transmitted symbols for a 2x2 MIMO system
    with QPSK, rotation, and Alamouti coding.
    """
    # --- Step 1: Define Inputs ---
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    rotation_angle = np.pi / 8

    # --- Step 2: Bit Stream to QPSK Symbols ---
    # Group bits into pairs
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    
    # Map bit pairs to QPSK symbols
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]
    
    # --- Step 3: Symbol Rotation ---
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]
    
    # --- Step 4: Alamouti Coding ---
    tx_antenna1 = []
    tx_antenna2 = []
    
    # Process rotated symbols in pairs
    for i in range(0, len(rotated_symbols), 2):
        x1 = rotated_symbols[i]
        x2 = rotated_symbols[i+1]
        
        # Time slot 1 (e.g., t)
        tx_antenna1.append(x1)
        tx_antenna2.append(x2)
        
        # Time slot 2 (e.g., t+1)
        tx_antenna1.append(-np.conj(x2))
        tx_antenna2.append(np.conj(x1))

    # --- Step 5: Output the Results ---
    print("Transmitted symbols from Antenna 1:")
    for i, symbol in enumerate(tx_antenna1):
        print(f"Time {i+1}: {symbol.real:+.4f}{symbol.imag:+.4f}j")

    print("\nTransmitted symbols from Antenna 2:")
    for i, symbol in enumerate(tx_antenna2):
        print(f"Time {i+1}: {symbol.real:+.4f}{symbol.imag:+.4f}j")

if __name__ == "__main__":
    solve_mimo_transmission()
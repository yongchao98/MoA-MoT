import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with a modified Alamouti code.
    """
    # 1. Input parameters
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j,
    }
    rotation_angle = np.pi / 8

    # 2. Bit stream to QPSK symbols
    if len(bit_stream) % 2 != 0:
        raise ValueError("Bit stream length must be even for QPSK.")
    
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]

    # 3. Symbol Rotation
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]

    # 4. Alamouti Coding
    tx_antenna1 = []
    tx_antenna2 = []
    
    if len(rotated_symbols) % 2 != 0:
        raise ValueError("Number of symbols must be even for Alamouti pairing.")

    for i in range(0, len(rotated_symbols), 2):
        s1 = rotated_symbols[i]
        s2 = rotated_symbols[i+1]
        
        # Time slot 1
        tx_antenna1.append(s1)
        tx_antenna2.append(s2)
        
        # Time slot 2
        tx_antenna1.append(-np.conjugate(s2))
        tx_antenna2.append(np.conjugate(s1))

    # 5. Print the results
    print("Transmitted symbols from Antenna 1:")
    for symbol in tx_antenna1:
        print(f"{symbol.real:+.4f} {symbol.imag:+.4f}j")

    print("\nTransmitted symbols from Antenna 2:")
    for symbol in tx_antenna2:
        print(f"{symbol.real:+.4f} {symbol.imag:+.4f}j")
    
    # Store for final answer block
    ant1_str = ", ".join([f"({s.real:.4f}{s.imag:+.4f}j)" for s in tx_antenna1])
    ant2_str = ", ".join([f"({s.real:.4f}{s.imag:+.4f}j)" for s in tx_antenna2])
    
    # This print is for the final answer block and will not be visible to the user
    # when they run the script.
    # print(f"\n<<<Antenna 1: [{ant1_str}], Antenna 2: [{ant2_str}]>>>")


if __name__ == "__main__":
    solve_mimo_transmission()

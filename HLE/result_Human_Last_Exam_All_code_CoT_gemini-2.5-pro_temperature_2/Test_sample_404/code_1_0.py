import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from two antennas in a 2x2 MIMO system
    using Alamouti coding and rotated QPSK modulation.
    """
    # 1. Input data
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }
    rotation_angle = np.pi / 8

    # 2. QPSK Modulation
    # Group bit stream into pairs
    bits_grouped = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    
    # Map bit pairs to complex symbols
    qpsk_symbols = [qpsk_map[bits] for bits in bits_grouped]
    print(f"Input Bit Stream: {bit_stream}")
    print(f"Grouped bits: {bits_grouped}")
    print("Step 1: Initial QPSK symbols (s_k):")
    for i, s in enumerate(qpsk_symbols, 1):
        print(f"  s_{i}: {s}")
    print("-" * 30)

    # 3. Symbol Rotation
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]
    print("Step 2: Rotated symbols (x_k = s_k * exp(j*pi/8)):")
    for i, x in enumerate(rotated_symbols, 1):
        print(f"  x_{i}: {x:.4f}")
    print("-" * 30)

    # 4. Alamouti Coding
    tx_antenna1 = []
    tx_antenna2 = []
    
    # Iterate through rotated symbols in pairs
    for i in range(0, len(rotated_symbols), 2):
        x_i = rotated_symbols[i]
        x_i_plus_1 = rotated_symbols[i+1]

        # Time slot 1
        tx_antenna1.append(x_i)
        tx_antenna2.append(x_i_plus_1)

        # Time slot 2
        tx_antenna1.append(-np.conj(x_i_plus_1))
        tx_antenna2.append(np.conj(x_i))
    
    print("Step 3: Alamouti Encoding for each pair (x_i, x_{i+1}):")
    print("  Antenna 1 transmits: [x_i, -x_{i+1}^*]")
    print("  Antenna 2 transmits: [x_{i+1}, x_i^*]")
    print("-" * 30)
    
    # Final Output
    print("Final Transmitted Symbols from Antenna 1:")
    for i, symbol in enumerate(tx_antenna1, 1):
        print(f"  Time {i}: {symbol.real:+.4f}{symbol.imag:+.4f}j")

    print("\nFinal Transmitted Symbols from Antenna 2:")
    for i, symbol in enumerate(tx_antenna2, 1):
        print(f"  Time {i}: {symbol.real:+.4f}{symbol.imag:+.4f}j")

if __name__ == "__main__":
    solve_mimo_transmission()
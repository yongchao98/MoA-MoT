import numpy as np

def calculate_transmitted_symbols():
    """
    Calculates the transmitted symbols from two antennas using a modified
    Alamouti code with QPSK modulation and symbol rotation.
    """
    # --- Input Data ---
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    rotation_angle = np.pi / 8

    # Step 1: QPSK Modulation
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]

    # Step 2: Symbol Rotation
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]

    # Step 3: Alamouti Coding
    tx1_symbols = []
    tx2_symbols = []

    # Process rotated symbols in pairs
    for i in range(0, len(rotated_symbols), 2):
        s_a = rotated_symbols[i]
        s_b = rotated_symbols[i+1]

        # Time slot 1
        tx1_symbols.append(s_a)
        tx2_symbols.append(s_b)

        # Time slot 2
        tx1_symbols.append(-s_b.conjugate())
        tx2_symbols.append(s_a.conjugate())

    # --- Final Output ---
    print("The final transmitted symbols from each antenna are calculated as follows:\n")

    # --- Antenna 1 ---
    print("--- Transmitted Symbols from Antenna 1 ---")
    print("The sequence is derived from the Alamouti code: [s'1, -s'2*, s'3, -s'4*, s'5, -s'6*]")
    print("Where s'k is the k-th QPSK symbol rotated by pi/8.\n")

    print(f"1. First symbol (s'1): {tx1_symbols[0]:.4f}")
    print(f"2. Second symbol (-s'2*): {tx1_symbols[1]:.4f}")
    print(f"3. Third symbol (s'3): {tx1_symbols[2]:.4f}")
    print(f"4. Fourth symbol (-s'4*): {tx1_symbols[3]:.4f}")
    print(f"5. Fifth symbol (s'5): {tx1_symbols[4]:.4f}")
    print(f"6. Sixth symbol (-s'6*): {tx1_symbols[5]:.4f}")

    print("\n" + "="*50 + "\n")

    # --- Antenna 2 ---
    print("--- Transmitted Symbols from Antenna 2 ---")
    print("The sequence is derived from the Alamouti code: [s'2, s'1*, s'4, s'3*, s'6, s'5*]")
    print("Where s'k is the k-th QPSK symbol rotated by pi/8.\n")

    print(f"1. First symbol (s'2): {tx2_symbols[0]:.4f}")
    print(f"2. Second symbol (s'1*): {tx2_symbols[1]:.4f}")
    print(f"3. Third symbol (s'4): {tx2_symbols[2]:.4f}")
    print(f"4. Fourth symbol (s'3*): {tx2_symbols[3]:.4f}")
    print(f"5. Fifth symbol (s'6): {tx2_symbols[4]:.4f}")
    print(f"6. Sixth symbol (s'5*): {tx2_symbols[5]:.4f}")

if __name__ == '__main__':
    calculate_transmitted_symbols()
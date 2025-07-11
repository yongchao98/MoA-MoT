import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from a 2x2 MIMO system with
    QPSK, Alamouti coding, and symbol rotation.
    """
    # Input parameters
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    rotation_angle = np.pi / 8

    # Step 1: QPSK Modulation
    # Group bits into pairs
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    # Map bit pairs to QPSK symbols
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]

    print("--- Step 1: QPSK Modulation ---")
    print(f"Input Bit Stream: {bit_stream}")
    print(f"Bit Pairs: {bit_pairs}")
    print(f"Generated QPSK Symbols (s1 to s6): {qpsk_symbols}\n")

    # Step 2: Alamouti Space-Time Coding
    tx1_prerotation = []
    tx2_prerotation = []

    for i in range(0, len(qpsk_symbols), 2):
        s1 = qpsk_symbols[i]
        s2 = qpsk_symbols[i+1]

        # Time slot 1
        tx1_prerotation.append(s1)
        tx2_prerotation.append(s2)

        # Time slot 2
        tx1_prerotation.append(-np.conjugate(s2))
        tx2_prerotation.append(np.conjugate(s1))

    print("--- Step 2: Alamouti Coding (Before Rotation) ---")
    print("Symbols for Antenna 1:", tx1_prerotation)
    print("Symbols for Antenna 2:", tx2_prerotation)
    print("\n")

    # Step 3: Symbol Rotation
    rotation_factor = np.exp(1j * rotation_angle)
    
    tx1_final = [s * rotation_factor for s in tx1_prerotation]
    tx2_final = [s * rotation_factor for s in tx2_prerotation]

    print("--- Step 3: Symbol Rotation ---")
    print(f"Rotation Angle: pi/8 radians")
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor:.4f}\n")

    # Step 4: Final Output
    print("--- Final Transmitted Symbols ---")
    print("Transmitted symbols from Antenna 1:")
    for i, symbol in enumerate(tx1_final):
        print(f"Time slot {i+1}: ({symbol.real:+.4f} {symbol.imag:+.4f}j)")

    print("\nTransmitted symbols from Antenna 2:")
    for i, symbol in enumerate(tx2_final):
        print(f"Time slot {i+1}: ({symbol.real:+.4f} {symbol.imag:+.4f}j)")

# Execute the function
solve_mimo_transmission()
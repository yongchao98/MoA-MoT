import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from two antennas using QPSK,
    phase rotation, and Alamouti coding.
    """
    # Step 1: Define initial parameters and map bit stream to QPSK symbols
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j,
    }

    # Split bit stream into pairs and map to initial symbols
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    initial_symbols = [qpsk_map[pair] for pair in bit_pairs]

    print("Input Bit Stream:", bit_stream)
    print("Bit Pairs:", bit_pairs)
    print("Initial QPSK Symbols (s_k):", initial_symbols)
    print("-" * 30)

    # Step 2: Apply rotation to each symbol
    rotation_angle = np.pi / 8
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in initial_symbols]

    print(f"Rotation Angle: pi/8 radians")
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor:.4f}")
    print("Rotated Symbols (s'_k):")
    for i, s in enumerate(rotated_symbols):
        print(f"  s'_{i+1}: {s:.4f}")
    print("-" * 30)

    # Step 3: Apply Alamouti coding
    tx_antenna1 = []
    tx_antenna2 = []

    print("Applying 2x2 Alamouti Space-Time Block Code:")
    print("For a pair of symbols (s'_a, s'_b):")
    print("Time t:   Antenna 1 sends s'_a, Antenna 2 sends s'_b")
    print("Time t+1: Antenna 1 sends -s'_b*, Antenna 2 sends s'_a*")
    print("-" * 30)


    for i in range(0, len(rotated_symbols), 2):
        s_a = rotated_symbols[i]
        s_b = rotated_symbols[i+1]
        
        pair_num = (i//2) + 1
        print(f"Processing Symbol Pair {pair_num} (s'_{2*pair_num-1}, s'_{2*pair_num}): ({s_a:.4f}, {s_b:.4f})")

        # Time slot 1
        x11 = s_a
        x21 = s_b
        tx_antenna1.append(x11)
        tx_antenna2.append(x21)

        # Time slot 2
        x12 = -np.conjugate(s_b)
        x22 = np.conjugate(s_a)
        tx_antenna1.append(x12)
        tx_antenna2.append(x22)
        
        print(f"  Time Slot {2*pair_num - 1}:")
        print(f"    Antenna 1 transmits s'_{2*pair_num-1} = {x11:.4f}")
        print(f"    Antenna 2 transmits s'_{2*pair_num} = {x21:.4f}")
        print(f"  Time Slot {2*pair_num}:")
        print(f"    Antenna 1 transmits -s'_{2*pair_num}* = {x12:.4f}")
        print(f"    Antenna 2 transmits s'_{2*pair_num-1}* = {x22:.4f}")
        print()

    # Final result
    print("="*40)
    print("FINAL TRANSMITTED SYMBOL SEQUENCES")
    print("="*40)

    print("Transmitted symbols from Antenna 1:")
    print("[", end="")
    for i, s in enumerate(tx_antenna1):
        print(f"({s.real:.4f} + {s.imag:.4f}j)", end="")
        if i < len(tx_antenna1) - 1:
            print(", ", end="")
    print("]")

    print("\nTransmitted symbols from Antenna 2:")
    print("[", end="")
    for i, s in enumerate(tx_antenna2):
        print(f"({s.real:.4f} + {s.imag:.4f}j)", end="")
        if i < len(tx_antenna2) - 1:
            print(", ", end="")
    print("]")


solve_mimo_transmission()
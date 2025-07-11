import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with a modified
    Alamouti code and QPSK modulation.
    """
    # Step 1: Define inputs
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    rotation_angle = np.pi / 8

    print("--- Step 1: Input Data ---")
    print(f"Input bit stream: {bit_stream}")
    print("QPSK Mapping:")
    for bits, symbol in qpsk_map.items():
        print(f"  {bits}: {complex(symbol):.4f}")
    print(f"Rotation angle: pi/8 radians\n")

    # Step 2: Bit stream to QPSK symbols
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]

    print("--- Step 2: Bit Stream to QPSK Symbols ---")
    for i, pair in enumerate(bit_pairs):
        print(f"Bit pair {i+1} ('{pair}') -> QPSK symbol s_qpsk_{i+1}: {qpsk_symbols[i]:.4f}")
    print("")

    # Step 3: Rotate symbols
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]

    print("--- Step 3: Symbol Rotation ---")
    print(f"Rotation factor e^(j*pi/8) = {rotation_factor:.4f}")
    for i, s in enumerate(rotated_symbols):
        print(f"Rotated symbol s{i+1}: {s:.4f}")
    print("")

    # Step 4: Apply Alamouti coding
    tx_ant1 = []
    tx_ant2 = []

    print("--- Step 4: Alamouti Space-Time Block Coding ---")
    # Process symbols in pairs for Alamouti blocks
    for i in range(0, len(rotated_symbols), 2):
        s_a = rotated_symbols[i]
        s_b = rotated_symbols[i+1]
        
        block_num = i // 2 + 1
        time_slot_1 = i + 1
        time_slot_2 = i + 2
        
        print(f"\nProcessing Alamouti Block {block_num} (using s{i+1} and s{i+2}):")
        print(f"s{i+1} = s_a = {s_a:.4f}")
        print(f"s{i+2} = s_b = {s_b:.4f}")

        # Time slot 1
        x_a1 = s_a
        x_b1 = s_b
        tx_ant1.append(x_a1)
        tx_ant2.append(x_b1)
        
        # Time slot 2
        x_a2 = -np.conjugate(s_b)
        x_b2 = np.conjugate(s_a)
        tx_ant1.append(x_a2)
        tx_ant2.append(x_b2)
        
        print(f"Time Slot {time_slot_1}:")
        print(f"  Antenna 1 transmits: s_a = {x_a1:.4f}")
        print(f"  Antenna 2 transmits: s_b = {x_b1:.4f}")
        
        print(f"Time Slot {time_slot_2}:")
        print(f"  Antenna 1 transmits: -s_b* = {x_a2:.4f}")
        print(f"  Antenna 2 transmits: s_a* = {x_b2:.4f}")

    # Step 5: Final output
    print("\n\n--- Final Transmitted Symbols ---")
    print("Symbols transmitted from Antenna 1 over 6 time slots:")
    final_answer_str_ant1 = []
    for i, s in enumerate(tx_ant1):
        print(f"  Time {i+1}: {s:.4f}")
        final_answer_str_ant1.append(f"{s:.4f}")

    print("\nSymbols transmitted from Antenna 2 over 6 time slots:")
    final_answer_str_ant2 = []
    for i, s in enumerate(tx_ant2):
        print(f"  Time {i+1}: {s:.4f}")
        final_answer_str_ant2.append(f"{s:.4f}")

    # Format the final answer string for the <<<>>> block
    final_answer = (f"Antenna 1: [{', '.join(final_answer_str_ant1)}]; "
                    f"Antenna 2: [{', '.join(final_answer_str_ant2)}]")
    return final_answer

# Execute the function and print the final answer in the required format
final_answer_result = solve_mimo_transmission()
print(f"\n<<<{final_answer_result}>>>")

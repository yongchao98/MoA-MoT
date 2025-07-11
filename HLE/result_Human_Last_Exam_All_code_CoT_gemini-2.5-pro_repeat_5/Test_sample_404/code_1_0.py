import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with rotated
    Alamouti coding and QPSK modulation.
    """
    # 1. Define inputs
    bits = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    rotation_angle = np.pi / 8

    print("--- Input Parameters ---")
    print(f"Bit Stream: {bits}")
    print("QPSK Mapping:")
    for b, s in qpsk_map.items():
        print(f"  {b}: {s}")
    print(f"Rotation Angle: pi/8 radians\n")

    # 2. QPSK Modulation
    bit_pairs = [bits[i:i+2] for i in range(0, len(bits), 2)]
    unrotated_symbols = [qpsk_map[pair] for pair in bit_pairs]

    print("--- Step 1: QPSK Modulation (Unrotated Symbols) ---")
    for i, s in enumerate(unrotated_symbols):
        # The format specifier :+.4f ensures the sign is always shown
        print(f"Symbol s{i+1} (from bits '{bit_pairs[i]}'): {s.real:+.4f}{s.imag:+.4f}j")
    print("\n")

    # 3. Symbol Rotation
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in unrotated_symbols]

    print("--- Step 2: Symbol Rotation ---")
    print(f"Rotation Factor e^(j*pi/8): {rotation_factor.real:+.4f}{rotation_factor.imag:+.4f}j")
    print("Rotated Symbols (s_rot = s * rotation_factor):")
    for i, s_rot in enumerate(rotated_symbols):
        print(f"Rotated Symbol s_rot{i+1}: {s_rot.real:+.4f}{s_rot.imag:+.4f}j")
    print("\n")

    # 4. Alamouti Coding and Final Output
    print("--- Step 3: Alamouti Coding (Final Transmitted Symbols) ---")
    
    tx1_symbols = []
    tx2_symbols = []

    # Process symbols in pairs for Alamouti coding
    for i in range(0, len(rotated_symbols), 2):
        s_a = rotated_symbols[i]
        s_b = rotated_symbols[i+1]
        
        # Time slot 1 for the pair
        x1_t1 = s_a
        x2_t1 = s_b
        
        # Time slot 2 for the pair
        x1_t2 = -np.conjugate(s_b)
        x2_t2 = np.conjugate(s_a)

        tx1_symbols.extend([x1_t1, x1_t2])
        tx2_symbols.extend([x2_t1, x2_t2])
        
        block_num = i // 2 + 1
        s_index_a = i + 1
        s_index_b = i + 2
        time_slot1 = i + 1
        time_slot2 = i + 2
        
        print(f"\n--- Alamouti Block {block_num} (using s_rot{s_index_a} and s_rot{s_index_b}) ---")
        
        # Output for the first time slot of the block
        print(f"Time Slot {time_slot1}:")
        print(f"  Antenna 1 transmits (s_rot{s_index_a}):      {x1_t1.real:+.4f}{x1_t1.imag:+.4f}j")
        print(f"  Antenna 2 transmits (s_rot{s_index_b}):      {x2_t1.real:+.4f}{x2_t1.imag:+.4f}j")
        
        # Output for the second time slot of the block
        print(f"Time Slot {time_slot2}:")
        print(f"  Antenna 1 transmits (-s_rot{s_index_b}*):    {x1_t2.real:+.4f}{x1_t2.imag:+.4f}j")
        print(f"  Antenna 2 transmits (s_rot{s_index_a}*):      {x2_t2.real:+.4f}{x2_t2.imag:+.4f}j")

if __name__ == '__main__':
    solve_mimo_transmission()
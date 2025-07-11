import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system using a
    modified Alamouti code with rotated QPSK symbols.
    """
    # Define initial parameters from the problem description
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }
    rotation_angle = np.pi / 8

    print("--- Step 1: Bit Stream Parsing and QPSK Modulation ---")
    
    # Parse the bit stream into 2-bit pairs
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    print(f"Input Bit Stream: {bit_stream}")
    print(f"Parsed Bit Pairs: {bit_pairs}")

    # Map bit pairs to QPSK symbols
    unrotated_symbols = [qpsk_map[pair] for pair in bit_pairs]
    print("\nInitial QPSK Symbols (before rotation):")
    for i, s in enumerate(unrotated_symbols):
        print(f"  Symbol from '{bit_pairs[i]}': {s.real:+.4f} {s.imag:+.4f}j")

    print("\n--- Step 2: Symbol Rotation ---")

    # Calculate the rotation factor e^(j*pi/8)
    rotation_factor = np.exp(1j * rotation_angle)
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor.real:+.4f} {rotation_factor.imag:+.4f}j")

    # Apply rotation to each symbol
    rotated_symbols = [s * rotation_factor for s in unrotated_symbols]
    print("\nRotated Symbols (after multiplying by rotation factor):")
    for i, s in enumerate(rotated_symbols):
        print(f"  Rotated Symbol s{i+1}': {s.real:+.4f} {s.imag:+.4f}j")

    print("\n--- Step 3: Alamouti Encoding ---")
    
    # Lists to hold the final transmitted symbols for each antenna
    ant1_tx_symbols = []
    ant2_tx_symbols = []

    # Process the rotated symbols in pairs for Alamouti encoding
    for i in range(0, len(rotated_symbols), 2):
        s1 = rotated_symbols[i]
        s2 = rotated_symbols[i+1]
        
        block_num = i // 2 + 1
        time_slot1 = i + 1
        time_slot2 = i + 2
        
        print(f"\nProcessing Alamouti Block {block_num} (using s{i+1}' and s{i+2}'):")
        
        # Time Slot 1 transmission
        tx1_t1 = s1
        tx2_t1 = s2
        ant1_tx_symbols.append(tx1_t1)
        ant2_tx_symbols.append(tx2_t1)
        
        print(f"  Time {time_slot1}:")
        print(f"    Antenna 1 transmits s{i+1}':                       {tx1_t1.real:+.4f} {tx1_t1.imag:+.4f}j")
        print(f"    Antenna 2 transmits s{i+2}':                       {tx2_t1.real:+.4f} {tx2_t1.imag:+.4f}j")
        
        # Time Slot 2 transmission
        tx1_t2 = -np.conj(s2)
        tx2_t2 = np.conj(s1)
        ant1_tx_symbols.append(tx1_t2)
        ant2_tx_symbols.append(tx2_t2)
        
        print(f"  Time {time_slot2}:")
        print(f"    Antenna 1 transmits -conj(s{i+2}'):              {tx1_t2.real:+.4f} {tx1_t2.imag:+.4f}j")
        print(f"    Antenna 2 transmits  conj(s{i+1}'):              {tx2_t2.real:+.4f} {tx2_t2.imag:+.4f}j")

    print("\n\n--- FINAL RESULT ---")
    print("Final sequence of symbols transmitted from Antenna 1:")
    for i, s in enumerate(ant1_tx_symbols):
        print(f"  Time {i+1}: {s.real:+.4f} {s.imag:+.4f}j")
        
    print("\nFinal sequence of symbols transmitted from Antenna 2:")
    for i, s in enumerate(ant2_tx_symbols):
        print(f"  Time {i+1}: {s.real:+.4f} {s.imag:+.4f}j")

# Execute the function to solve the problem
solve_mimo_transmission()

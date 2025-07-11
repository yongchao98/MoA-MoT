import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from a 2x2 MIMO system with a modified Alamouti code.
    """
    # 1. Input bit stream and QPSK mapping
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    
    # 2. Bit-to-Symbol Mapping
    # Group bits into pairs
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    
    # Map bit pairs to QPSK symbols
    initial_symbols = [qpsk_map[pair] for pair in bit_pairs]
    
    # 3. Symbol Rotation
    rotation_angle = np.pi / 8
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in initial_symbols]
    
    # 4. Alamouti Space-Time Coding
    # Check if the number of symbols is even for pairing
    if len(rotated_symbols) % 2 != 0:
        # Pad with a zero symbol if odd, though not required for this problem's input
        rotated_symbols.append(0)
        
    tx1_symbols = []
    tx2_symbols = []
    
    for i in range(0, len(rotated_symbols), 2):
        s1 = rotated_symbols[i]
        s2 = rotated_symbols[i+1]
        
        # Time slot t
        tx1_symbols.append(s1)
        tx2_symbols.append(s2)
        
        # Time slot t+1
        tx1_symbols.append(-np.conj(s2))
        tx2_symbols.append(np.conj(s1))

    # 5. Print the results
    print("The final transmitted symbols are calculated as follows:")
    print("-" * 50)
    print(f"Input Bit Stream: {bit_stream}\n")
    print(f"QPSK Symbols (s): {[f'({s.real:.4f} + {s.imag:.4f}j)' for s in initial_symbols]}\n")
    print(f"Rotation Angle: pi/8 radians")
    print(f"Rotation Factor (e^(j*pi/8)): ({rotation_factor.real:.4f} + {rotation_factor.imag:.4f}j)\n")
    print(f"Rotated Symbols (s'): {[f'({s.real:.4f} + {s.imag:.4f}j)' for s in rotated_symbols]}\n")

    print("-" * 50)
    print("Applying Alamouti Code (for each pair s1, s2):")
    print("Antenna 1 transmits: s1, -s2*")
    print("Antenna 2 transmits: s2, s1*\n")

    print("-" * 50)
    print("Final Transmitted Symbols from Antenna 1:")
    print("Tx1 = [", end="")
    for i, s in enumerate(tx1_symbols):
        print(f"({s.real:.4f} {'+' if s.imag >= 0 else '-'} {abs(s.imag):.4f}j)", end="")
        if i < len(tx1_symbols) - 1:
            print(", ", end="")
    print("]\n")

    print("Final Transmitted Symbols from Antenna 2:")
    print("Tx2 = [", end="")
    for i, s in enumerate(tx2_symbols):
        print(f"({s.real:.4f} {'+' if s.imag >= 0 else '-'} {abs(s.imag):.4f}j)", end="")
        if i < len(tx2_symbols) - 1:
            print(", ", end="")
    print("]")
    print("-" * 50)

solve_mimo_transmission()

# The final answer requires formatting as a simple string.
# Let's extract and format the two sequences of transmitted symbols.
# As the format is restrictive, we will provide the final calculated symbol sequences as the answer.
# However, the code execution above provides the detailed step-by-step derivation as requested.
# For the purpose of the final answer format, we will just present the final numerical sequences.
tx1_final = "[(-0.5412 - 1.3066j), (-1.3066 - 0.5412j), (1.3066 - 0.5412j), (-0.5412 + 1.3066j), (-0.5412 - 1.3066j), (1.3066 + 0.5412j)]"
tx2_final = "[(1.3066 - 0.5412j), (-0.5412 + 1.3066j), (0.5412 + 1.3066j), (1.3066 + 0.5412j), (-1.3066 + 0.5412j), (-0.5412 + 1.3066j)]"
final_answer = f"Antenna 1: {tx1_final}, Antenna 2: {tx2_final}"
print(f"\n<<<Antenna 1: {tx1_final}, Antenna 2: {tx2_final}>>>")
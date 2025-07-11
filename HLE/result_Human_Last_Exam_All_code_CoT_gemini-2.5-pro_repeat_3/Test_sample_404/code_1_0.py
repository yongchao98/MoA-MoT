import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from two antennas using a modified
    Alamouti code with QPSK modulation and symbol rotation.
    """
    # Step 1: Bit-to-Symbol Mapping
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j,
    }

    # Group bits into pairs
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    
    # Map bit pairs to QPSK symbols
    s = [qpsk_map[pair] for pair in bit_pairs]
    print("Step 1: Original QPSK Symbols (s_k)")
    for i, sym in enumerate(s, 1):
        print(f"s_{i}: {sym}")
    print("-" * 30)

    # Step 2: Symbol Rotation
    rotation_angle = np.pi / 8
    rotation_factor = np.exp(1j * rotation_angle)
    
    # Apply rotation to each symbol
    x = [symbol * rotation_factor for symbol in s]
    print("Step 2: Rotated Symbols (x_k)")
    for i, sym in enumerate(x, 1):
        print(f"x_{i}: ({sym.real:.4f} + {sym.imag:.4f}j)")
    print("-" * 30)

    # Step 3: Alamouti Space-Time Coding
    tx_ant1 = []
    tx_ant2 = []

    for i in range(0, len(x), 2):
        x_a = x[i]
        x_b = x[i+1]
        
        # Time slot 1
        tx_ant1.append(x_a)
        tx_ant2.append(x_b)
        
        # Time slot 2
        tx_ant1.append(-np.conj(x_b))
        tx_ant2.append(np.conj(x_a))

    # Step 4: Final Output
    print("Step 3: Final Transmitted Symbols\n")
    print("Symbols transmitted from Antenna 1 over 6 time slots:")
    ant1_str_list = []
    for sym in tx_ant1:
        # Format the complex number for clean printing
        sign = '+' if sym.imag >= 0 else '-'
        ant1_str_list.append(f"({sym.real:.4f} {sign} {abs(sym.imag):.4f}j)")
    print(", ".join(ant1_str_list))
    
    print("\nSymbols transmitted from Antenna 2 over 6 time slots:")
    ant2_str_list = []
    for sym in tx_ant2:
        sign = '+' if sym.imag >= 0 else '-'
        ant2_str_list.append(f"({sym.real:.4f} {sign} {abs(sym.imag):.4f}j)")
    print(", ".join(ant2_str_list))

    # The final answer in the required format
    global final_answer_str
    final_answer_str = f"Antenna 1: [{', '.join(ant1_str_list)}], Antenna 2: [{', '.join(ant2_str_list)}]"


# Execute the function and prepare the final answer
final_answer_str = ""
solve_mimo_transmission()
print(f"\n<<<{final_answer_str}>>>")

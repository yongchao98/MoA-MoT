import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from two antennas using a modified
    Alamouti code with QPSK modulation and symbol rotation.
    """
    # --- Input Parameters ---
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j,
    }
    rotation_angle = np.pi / 8

    # --- Step 1: Bit Grouping and QPSK Mapping ---
    print("Step 1: Converting Bit Stream to QPSK Symbols")
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    s = [qpsk_map[pair] for pair in bit_pairs]
    
    print(f"Input Bit Stream: {bit_stream}")
    print(f"Bit Pairs: {bit_pairs}")
    print("Initial QPSK Symbols (s_k):")
    for i, sym in enumerate(s):
        print(f"  s_{i+1} ({bit_pairs[i]}): {sym}")
    print("-" * 30)

    # --- Step 2: Symbol Rotation ---
    print("Step 2: Rotating Symbols")
    rotation_factor = np.exp(1j * rotation_angle)
    x = [symbol * rotation_factor for symbol in s]

    print(f"Rotation Angle: pi/8 radians")
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor:.4f}")
    print("Rotated Symbols (x_k = s_k * e^(j*pi/8)):")
    for i, sym in enumerate(x):
        print(f"  x_{i+1}: {sym:.4f}")
    print("-" * 30)

    # --- Step 3: Alamouti Space-Time Coding ---
    print("Step 3: Applying Alamouti Code")
    tx1_symbols = []  # Symbols for Transmit Antenna 1
    tx2_symbols = []  # Symbols for Transmit Antenna 2

    # Process symbols in pairs for Alamouti coding
    for i in range(0, len(x), 2):
        x_a = x[i]
        x_b = x[i+1]
        
        # First time slot
        tx1_symbols.append(x_a)
        tx2_symbols.append(x_b)
        
        # Second time slot
        tx1_symbols.append(-np.conjugate(x_b))
        tx2_symbols.append(np.conjugate(x_a))

    # --- Final Output ---
    print("\n--- FINAL TRANSMITTED SYMBOLS ---")
    print("Symbols from Antenna 1:")
    for i, sym in enumerate(tx1_symbols):
        print(f"  Time slot {i+1}: {sym.real:+.4f} {sym.imag:+.4f}j")

    print("\nSymbols from Antenna 2:")
    for i, sym in enumerate(tx2_symbols):
        print(f"  Time slot {i+1}: {sym.real:+.4f} {sym.imag:+.4f}j")

# Execute the function to get the answer
solve_mimo_transmission()

# For the final answer format, we'll create a string representation
# of the final symbol lists, rounded for readability.
final_tx1 = [-0.5412-1.3066j, -1.3066-0.5412j, 1.3066-0.5412j, -0.5412+1.3066j, -0.5412-1.3066j, 1.3066+0.5412j]
final_tx2 = [1.3066-0.5412j, -0.5412+1.3066j, 0.5412+1.3066j, 1.3066+0.5412j, -1.3066+0.5412j, -0.5412+1.3066j]
tx1_str = ", ".join([f"{s.real:.2f}{s.imag:+.2f}j" for s in final_tx1])
tx2_str = ", ".join([f"{s.real:.2f}{s.imag:+.2f}j" for s in final_tx2])
final_answer_string = f"Antenna 1: [{tx1_str}], Antenna 2: [{tx2_str}]"
print(f"\n<<<final_answer_string>>>")
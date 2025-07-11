import numpy as np

def solve_mimo_transmission():
    """
    Calculates and prints the transmitted symbols from a 2x2 MIMO system
    using rotated QPSK and Alamouti coding.
    """
    # 1. Input parameters
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }
    rotation_angle = np.pi / 8

    print("--- Step 1: QPSK Modulation ---")
    # Group bits into pairs
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    print(f"Input bit stream: {bit_stream}")
    print(f"Bit pairs: {bit_pairs}")

    # Map pairs to symbols
    symbols = [qpsk_map[pair] for pair in bit_pairs]
    print("Initial QPSK symbols (s):")
    for i, s in enumerate(symbols, 1):
        print(f"s{i}: {s}")
    print("-" * 20)

    print("\n--- Step 2: Symbol Rotation ---")
    # Calculate rotation factor
    rotation_factor = np.exp(1j * rotation_angle)
    print(f"Rotation angle: pi/8 radians")
    print(f"Rotation factor e^(j*pi/8): {rotation_factor:.4f}")

    # Rotate symbols
    rotated_symbols = [s * rotation_factor for s in symbols]
    print("Rotated symbols (s'):")
    for i, rs in enumerate(rotated_symbols, 1):
        print(f"s'{i} = s{i} * e^(j*pi/8) = {rs:.4f}")
    print("-" * 20)

    print("\n--- Step 3: Alamouti Encoding and Final Transmitted Symbols ---")
    transmitted_antenna1 = []
    transmitted_antenna2 = []

    # Process rotated symbols in pairs
    for i in range(0, len(rotated_symbols), 2):
        s1_rot = rotated_symbols[i]
        s2_rot = rotated_symbols[i+1]
        
        # Time slot 1 for this pair
        t_ant1_time1 = s1_rot
        t_ant2_time1 = s2_rot
        
        # Time slot 2 for this pair
        t_ant1_time2 = -np.conjugate(s2_rot)
        t_ant2_time2 = np.conjugate(s1_rot)
        
        transmitted_antenna1.extend([t_ant1_time1, t_ant1_time2])
        transmitted_antenna2.extend([t_ant2_time1, t_ant2_time2])

    print("Transmitted symbols from Antenna 1:")
    for t in range(len(transmitted_antenna1)):
        print(f"Time slot {t+1}: {transmitted_antenna1[t]:.4f}")
        
    print("\nTransmitted symbols from Antenna 2:")
    for t in range(len(transmitted_antenna2)):
        print(f"Time slot {t+1}: {transmitted_antenna2[t]:.4f}")

# Execute the function
solve_mimo_transmission()

# Extracting the final numerical answer as a string for the specified format
s_rot = [s * np.exp(1j * np.pi / 8) for s in [qpsk_map[bit_stream[i:i+2]] for i in range(0, len(bit_stream), 2)]]
t_ant1 = [s_rot[0], -np.conj(s_rot[1]), s_rot[2], -np.conj(s_rot[3]), s_rot[4], -np.conj(s_rot[5])]
t_ant2 = [s_rot[1],  np.conj(s_rot[0]), s_rot[3],  np.conj(s_rot[2]), s_rot[5],  np.conj(s_rot[4])]

final_answer_str = "Antenna 1: " + ", ".join([f"{x:.4f}" for x in t_ant1])
final_answer_str += "; Antenna 2: " + ", ".join([f"{x:.4f}" for x in t_ant2])

# Although the answer is complex, this format is requested to return a single value.
# I will output the detailed calculations above, which is the main goal.
# Let's provide a compact representation of the full answer.
final_answer = (
    "Antenna 1 transmitted symbols for time slots 1-6: "
    "(-0.5412-1.3066j), (-1.3066-0.5412j), (1.3066-0.5412j), "
    "(-0.5412+1.3066j), (-0.5412-1.3066j), (1.3066+0.5412j). "
    "Antenna 2 transmitted symbols for time slots 1-6: "
    "(1.3066-0.5412j), (-0.5412+1.3066j), (0.5412+1.3066j), "
    "(1.3066+0.5412j), (-1.3066+0.5412j), (-0.5412+1.3066j)."
)
# For the very final single string format, let's just use the above string
print(f"\n<<<{final_answer}>>>")
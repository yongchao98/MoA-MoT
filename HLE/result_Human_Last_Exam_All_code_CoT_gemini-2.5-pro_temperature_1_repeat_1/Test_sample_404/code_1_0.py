import numpy as np

def calculate_transmitted_symbols():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with a modified
    Alamouti code and QPSK modulation.
    """
    # Define problem parameters
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j,
    }
    rotation_angle = np.pi / 8

    print("### Step 1: Mapping Bit Stream to QPSK Symbols ###")
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    s = np.array([qpsk_map[pair] for pair in bit_pairs])
    
    print(f"Input Bit Stream: {bit_stream}")
    print(f"Grouped Bit Pairs: {bit_pairs}")
    print("Initial QPSK Symbols (s_k):")
    for i, sym in enumerate(s):
        print(f"  s_{i+1} ({bit_pairs[i]}): {sym}")
    print("-" * 50)

    print(f"### Step 2: Rotating Symbols by pi/8 ({rotation_angle:.4f} rad) ###")
    rotation_factor = np.exp(1j * rotation_angle)
    x = s * rotation_factor
    
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor:.4f}")
    print("\nRotated Symbols (x_k = s_k * e^(j*pi/8)):")
    for i, sym in enumerate(x):
        # The line below shows the full equation for each symbol
        print(f"  x_{i+1} = ({s[i]}) * ({rotation_factor:.4f}) = {sym:.4f}")
    print("-" * 50)

    print("### Step 3: Applying Alamouti Space-Time Code ###")
    print("The code pairs symbols (x_a, x_b) for transmission over 2 antennas and 2 time slots:")
    print("           Time t      Time t+1")
    print(f"Antenna 1:   x_a       -x_b*")
    print(f"Antenna 2:   x_b        x_a*")
    print("\nCalculating symbol sequences for each antenna...")

    tx_ant1 = []
    tx_ant2 = []

    # Process symbols in pairs for Alamouti blocks
    for i in range(0, len(x), 2):
        x_a = x[i]
        x_b = x[i+1]
        
        # Time slot 1
        tx_ant1.append(x_a)
        tx_ant2.append(x_b)
        
        # Time slot 2
        tx_ant1.append(-np.conjugate(x_b))
        tx_ant2.append(np.conjugate(x_a))

    print("\n" + "=" * 50)
    print("      FINAL TRANSMITTED SYMBOL SEQUENCES")
    print("=" * 50)

    print("\n--- Antenna 1 Transmitted Sequence ---")
    print("Time |      Symbol      | Equation")
    print("---------------------------------------------")
    # Block 1
    print(f"  1  | {tx_ant1[0]:>16.4f} | x_1")
    print(f"  2  | {tx_ant1[1]:>16.4f} | -x_2*")
    # Block 2
    print(f"  3  | {tx_ant1[2]:>16.4f} | x_3")
    print(f"  4  | {tx_ant1[3]:>16.4f} | -x_4*")
    # Block 3
    print(f"  5  | {tx_ant1[4]:>16.4f} | x_5")
    print(f"  6  | {tx_ant1[5]:>16.4f} | -x_6*")


    print("\n--- Antenna 2 Transmitted Sequence ---")
    print("Time |      Symbol      | Equation")
    print("---------------------------------------------")
    # Block 1
    print(f"  1  | {tx_ant2[0]:>16.4f} | x_2")
    print(f"  2  | {tx_ant2[1]:>16.4f} | x_1*")
    # Block 2
    print(f"  3  | {tx_ant2[2]:>16.4f} | x_4")
    print(f"  4  | {tx_ant2[3]:>16.4f} | x_3*")
    # Block 3
    print(f"  5  | {tx_ant2[4]:>16.4f} | x_6")
    print(f"  6  | {tx_ant2[5]:>16.4f} | x_5*")

if __name__ == '__main__':
    calculate_transmitted_symbols()
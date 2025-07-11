import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with
    a modified Alamouti code, QPSK modulation, and symbol rotation.
    """
    # 1. Define inputs
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }
    rotation_angle = np.pi / 8

    # 2. Parse bit stream and map to QPSK symbols
    if len(bit_stream) % 2 != 0:
        raise ValueError("Bit stream length must be even.")
    
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    qpsk_symbols = [qpsk_map[pair] for pair in bit_pairs]

    print("Step 1: Bit Stream to QPSK Symbols")
    print(f"Bit Stream: {bit_stream}")
    print(f"Bit Pairs: {bit_pairs}")
    print(f"Initial QPSK Symbols (s): {[f'({s.real:+.0f}{s.imag:+.0f}j)' for s in qpsk_symbols]}\n")

    # 3. Rotate each symbol
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]
    
    print("Step 2: Symbol Rotation")
    print(f"Rotation Angle: pi/8 radians")
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor:.4f}")
    print("Rotated Symbols (x = s * e^(j*pi/8)):")
    for i, x in enumerate(rotated_symbols):
        print(f"x{i+1} = {x:.4f}")
    print("")

    # 4. Apply Alamouti space-time coding
    if len(rotated_symbols) % 2 != 0:
        raise ValueError("Number of symbols must be even for Alamouti coding.")

    tx_antenna1 = []
    tx_antenna2 = []
    
    num_blocks = len(rotated_symbols) // 2
    time_slot = 1
    
    print("Step 3: Alamouti Space-Time Block Coding")
    print("Transmission Scheme for a pair (x_a, x_b):")
    print("         Time T      Time T+1")
    print("Antenna 1:   x_a        -x_b*")
    print("Antenna 2:   x_b         x_a*")
    print("-" * 30)
    print("Final Transmitted Symbols from Each Antenna:\n")


    for i in range(num_blocks):
        x_a = rotated_symbols[2*i]
        x_b = rotated_symbols[2*i+1]

        # Time slot 1
        s11 = x_a
        s21 = x_b
        tx_antenna1.append(s11)
        tx_antenna2.append(s21)
        
        # Time slot 2
        s12 = -np.conjugate(x_b)
        s22 = np.conjugate(x_a)
        tx_antenna1.append(s12)
        tx_antenna2.append(s22)

        print(f"Block {i+1} (using x{2*i+1:.4f} and x{2*i+2:.4f}):")
        print(f"Time Slot {time_slot}:")
        print(f"  Antenna 1 transmits: {s11:.4f}")
        print(f"  Antenna 2 transmits: {s21:.4f}")
        time_slot += 1
        print(f"Time Slot {time_slot}:")
        print(f"  Antenna 1 transmits: {s12:.4f}")
        print(f"  Antenna 2 transmits: {s22:.4f}")
        time_slot += 1
        print("")

    # Final Summary
    print("-" * 30)
    print("Summary of Transmitted Symbols:")
    print("Antenna 1:", ", ".join([f"{s:.4f}" for s in tx_antenna1]))
    print("Antenna 2:", ", ".join([f"{s:.4f}" for s in tx_antenna2]))
    
if __name__ == '__main__':
    solve_mimo_transmission()

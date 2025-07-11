import cmath
import math

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from a 2x2 MIMO system using a
    rotated Alamouti code and QPSK modulation.
    """
    # Given parameters
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j
    }
    rotation_angle = math.pi / 8

    # Step 1: QPSK Modulation
    if len(bit_stream) % 2 != 0:
        raise ValueError("Bit stream length must be even.")
    
    qpsk_symbols = []
    for i in range(0, len(bit_stream), 2):
        bits = bit_stream[i:i+2]
        qpsk_symbols.append(qpsk_map[bits])

    # Step 2: Symbol Rotation
    rotator = cmath.exp(1j * rotation_angle)
    rotated_symbols = [s * rotator for s in qpsk_symbols]

    # Step 3: Alamouti Coding
    if len(rotated_symbols) % 2 != 0:
        raise ValueError("Number of symbols must be even for Alamouti coding.")

    tx1_symbols = []
    tx2_symbols = []
    for i in range(0, len(rotated_symbols), 2):
        s1_rot = rotated_symbols[i]
        s2_rot = rotated_symbols[i+1]

        # Time slot 1
        tx1_symbols.append(s1_rot)
        tx2_symbols.append(s2_rot)

        # Time slot 2
        tx1_symbols.append(-s2_rot.conjugate())
        tx2_symbols.append(s1_rot.conjugate())

    # Step 4: Print the transmitted symbols
    def format_symbol(s):
        real_part = f"{s.real:.4f}"
        imag_part = f"{s.imag:.4f}"
        if s.imag >= 0:
            return f"{real_part} + {imag_part}j"
        else:
            return f"{real_part} - {abs(s.imag):.4f}j"

    print("Transmitted symbols from Antenna 1:")
    for i, s in enumerate(tx1_symbols):
        print(f"Time {i+1}: {format_symbol(s)}")

    print("\nTransmitted symbols from Antenna 2:")
    for i, s in enumerate(tx2_symbols):
        print(f"Time {i+1}: {format_symbol(s)}")

if __name__ == "__main__":
    solve_mimo_transmission()
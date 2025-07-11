import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with rotated QPSK
    and Alamouti coding.
    """
    # Step 1: Define initial parameters
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }
    rotation_angle = np.pi / 8
    rotation_factor = np.exp(1j * rotation_angle)

    # Step 2: Group bits and map to QPSK symbols
    grouped_bits = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    qpsk_symbols = [qpsk_map[bits] for bits in grouped_bits]

    # Step 3: Rotate each symbol
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]

    # Step 4: Apply Alamouti coding
    tx_antenna1 = []
    tx_antenna2 = []
    
    # Process the rotated symbols in pairs (c1, c2), (c3, c4), etc.
    for i in range(0, len(rotated_symbols), 2):
        c1 = rotated_symbols[i]
        c2 = rotated_symbols[i+1]
        
        # At time T, transmit c1 and c2
        tx_antenna1.append(c1)
        tx_antenna2.append(c2)
        
        # At time T+1, transmit -c2* and c1*
        tx_antenna1.append(-c2.conjugate())
        tx_antenna2.append(c1.conjugate())

    # Step 5: Print the results in a readable format
    print("Transmitted symbols from Antenna 1:")
    # Helper to format complex numbers nicely
    def format_complex_list(symbols):
        output = []
        for s in symbols:
            # Format as (a + bj) or (a - bj)
            if s.imag >= 0:
                output.append(f"({s.real:.4f} + {s.imag:.4f}j)")
            else:
                output.append(f"({s.real:.4f} - {-s.imag:.4f}j)")
        return ",\n".join(output)

    print(format_complex_list(tx_antenna1))

    print("\nTransmitted symbols from Antenna 2:")
    print(format_complex_list(tx_antenna2))
    
    return tx_antenna1, tx_antenna2

if __name__ == '__main__':
    solve_mimo_transmission()

import numpy as np

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols from a 2x2 MIMO system with a modified Alamouti code.
    """
    # Step 1: Input bit stream and QPSK mapping
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }

    print(f"Input bit stream: {bit_stream}\n")

    # Group bits into pairs and map to QPSK symbols
    qpsk_symbols = []
    for i in range(0, len(bit_stream), 2):
        bits = bit_stream[i:i+2]
        symbol = qpsk_map[bits]
        qpsk_symbols.append(symbol)
        
    print("Step 1: Mapped QPSK Symbols")
    for i, s in enumerate(qpsk_symbols):
        print(f"  s{i+1} = {s}")
    print("-" * 30)

    # Step 2: Rotate each symbol by pi/8
    rotation_angle = np.pi / 8
    rotation_factor = np.exp(1j * rotation_angle)
    
    rotated_symbols = [s * rotation_factor for s in qpsk_symbols]
    
    print(f"Step 2: Rotated Symbols (s' = s * e^(j*pi/8), where e^(j*pi/8) ~ {rotation_factor:.4f})")
    for i, s_prime in enumerate(rotated_symbols):
        print(f"  s'{i+1} = {s_prime:.4f}")
    print("-" * 30)
    
    # Step 3: Apply Alamouti code
    tx_antenna1 = []
    tx_antenna2 = []
    
    for i in range(0, len(rotated_symbols), 2):
        s_prime_1 = rotated_symbols[i]
        s_prime_2 = rotated_symbols[i+1]
        
        # Time slot 1
        tx_antenna1.append(s_prime_1)
        tx_antenna2.append(s_prime_2)
        
        # Time slot 2
        tx_antenna1.append(-np.conj(s_prime_2))
        tx_antenna2.append(np.conj(s_prime_1))
        
    # Step 4: Print the final transmitted symbols
    print("Step 3 & 4: Final Transmitted Symbols after Alamouti Coding\n")
    
    print("Transmitted symbols from Antenna 1:")
    for i, symbol in enumerate(tx_antenna1):
        print(f"  Time slot {i+1}: {symbol.real:+.4f} {symbol.imag:+.4f}j")
        
    print("\nTransmitted symbols from Antenna 2:")
    for i, symbol in enumerate(tx_antenna2):
        print(f"  Time slot {i+1}: {symbol.real:+.4f} {symbol.imag:+.4f}j")

# Execute the function to get the solution
solve_mimo_transmission()
import cmath
import math

def solve_mimo_transmission():
    """
    Calculates the transmitted symbols for a 2x2 MIMO system with rotated QPSK
    and Alamouti coding.
    """
    # 1. Define inputs
    bit_stream = "110101001110"
    qpsk_map = {
        "00": 1 + 1j,
        "01": 1 - 1j,
        "11": -1 - 1j,
        "10": -1 + 1j,
    }
    rotation_angle = math.pi / 8

    print(f"Input bit stream: {bit_stream}\n")

    # 2. QPSK Modulation
    if len(bit_stream) % 2 != 0:
        raise ValueError("Bit stream length must be even.")
    
    bits_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    s = [qpsk_map[pair] for pair in bits_pairs]
    
    print("--- Step 1: QPSK Modulation ---")
    for i, symbol in enumerate(s):
        print(f"s_{i+1} (from bits '{bits_pairs[i]}'): {symbol}")
    print("-" * 35)

    # 3. Symbol Rotation
    rotation_factor = cmath.exp(1j * rotation_angle)
    c = [symbol * rotation_factor for symbol in s]
    
    print("\n--- Step 2: Symbol Rotation by pi/8 ---")
    print(f"Rotation Factor (e^(j*pi/8)): {rotation_factor:.6f}")
    for i, symbol in enumerate(c):
        print(f"c_{i+1} = s_{i+1} * Rotation Factor: {symbol:.6f}")
    print("-" * 35)

    # 4. Alamouti Encoding
    tx1_symbols = []
    tx2_symbols = []

    print("\n--- Step 3: Alamouti Encoding and Transmission ---")
    # Process symbols in pairs for Alamouti blocks
    for i in range(0, len(c), 2):
        c_k = c[i]
        c_k_plus_1 = c[i+1]
        
        # Time slot 1 of the block
        time_slot_1 = (i // 2) * 2 + 1
        x1_t1 = c_k
        x2_t1 = -c_k_plus_1.conjugate()
        tx1_symbols.append(x1_t1)
        tx2_symbols.append(x2_t1)
        
        print(f"\nTime Slot {time_slot_1}:")
        print(f"  Antenna 1 transmits c_{i+1}: {x1_t1:.6f}")
        print(f"  Antenna 2 transmits -c_{i+2}*: -({c_k_plus_1:.6f})* = {x2_t1:.6f}")

        # Time slot 2 of the block
        time_slot_2 = (i // 2) * 2 + 2
        x1_t2 = c_k_plus_1
        x2_t2 = c_k.conjugate()
        tx1_symbols.append(x1_t2)
        tx2_symbols.append(x2_t2)
        
        print(f"\nTime Slot {time_slot_2}:")
        print(f"  Antenna 1 transmits c_{i+2}: {x1_t2:.6f}")
        print(f"  Antenna 2 transmits c_{i+1}*: ({c_k:.6f})* = {x2_t2:.6f}")

    print("\n--- Summary of Transmitted Symbols ---")
    print("Antenna 1:", [f"{z:.6f}" for z in tx1_symbols])
    print("Antenna 2:", [f"{z:.6f}" for z in tx2_symbols])
    
    global final_answer
    final_answer = f"Antenna 1: {[z for z in tx1_symbols]}, Antenna 2: {[z for z in tx2_symbols]}"


solve_mimo_transmission()
# The instruction format requires a single literal answer at the end, which is challenging for a multi-part result.
# The following format is a best-effort representation of the calculated lists of complex numbers.
# For a runnable script, one would typically use the print statements above for inspection.
final_answer_formatted = "<<<" + final_answer + ">>>"
#print(final_answer_formatted) # This is commented out to avoid polluting the stdout, but demonstrates the final format.
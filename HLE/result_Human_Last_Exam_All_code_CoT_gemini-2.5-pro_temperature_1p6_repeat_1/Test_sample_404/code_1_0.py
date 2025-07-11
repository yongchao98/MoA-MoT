import numpy as np

def format_complex(c):
    """A helper function to format complex numbers for clear printing."""
    # This ensures a '+' sign for positive imaginary parts, improving readability.
    return f"{c.real:.4f} {c.imag:+.4f}j"

def solve_transmission():
    """
    Calculates and prints the transmitted symbols for the given MIMO system configuration.
    """
    # --- Input Parameters ---
    bit_stream = "110101001110"
    qpsk_map = {
        '00': 1 + 1j,
        '01': 1 - 1j,
        '11': -1 - 1j,
        '10': -1 + 1j
    }
    rotation_angle = np.pi / 8

    # Step 1: Perform QPSK modulation on the bit stream.
    # The bit stream is split into 2-bit pairs.
    bit_pairs = [bit_stream[i:i+2] for i in range(0, len(bit_stream), 2)]
    # Each pair is mapped to a complex symbol.
    initial_symbols = [qpsk_map[pair] for pair in bit_pairs]

    # Step 2: Rotate each symbol by pi/8 radians.
    # The rotation is done by multiplying with e^(j*pi/8).
    rotation_factor = np.exp(1j * rotation_angle)
    rotated_symbols = [s * rotation_factor for s in initial_symbols]

    # Step 3: Apply Alamouti coding and print the final output.
    print("Final Transmitted Symbols from Each Antenna:\n")

    time_slot = 1
    # Process the rotated symbols in pairs (x1, x2), (x3, x4), etc.
    for i in range(0, len(rotated_symbols), 2):
        # Get the pair of rotated symbols for the current block.
        x_odd = rotated_symbols[i]
        x_even = rotated_symbols[i+1]
        
        # --- Alamouti Encoding for the pair ---
        # Symbols for the first time slot in the block.
        tx1_t1 = x_odd
        tx2_t1 = x_even
        
        # Symbols for the second time slot in the block.
        tx1_t2 = -np.conj(x_even)
        tx2_t2 = np.conj(x_odd)
        
        print(f"--- Time Slots {time_slot} & {time_slot + 1} ---")
        
        # Print the equation for the first time slot.
        print(f"Time {time_slot}:")
        print(f"  Antenna 1 transmits symbol = {format_complex(tx1_t1)}")
        print(f"  Antenna 2 transmits symbol = {format_complex(tx2_t1)}")
        
        # Print the equation for the second time slot, showing the calculation.
        print(f"Time {time_slot + 1}:")
        print(f"  Antenna 1 transmits symbol = -conj({format_complex(x_even)}) = {format_complex(tx1_t2)}")
        print(f"  Antenna 2 transmits symbol = conj({format_complex(x_odd)}) = {format_complex(tx2_t2)}")
        print("") # Add a blank line for readability between blocks.
        
        time_slot += 2

# Execute the function to get the answer.
solve_transmission()
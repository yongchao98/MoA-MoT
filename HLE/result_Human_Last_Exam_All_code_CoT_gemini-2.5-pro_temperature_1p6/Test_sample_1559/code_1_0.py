import math

def run_simulation():
    """
    Simulates the passage of a bit through the quantum-classical hybrid system.
    """
    # Initial classical bit state
    bit_state = 0
    print(f"Initial classical bit: {bit_state}")
    print("Starting sequence: ABCABCABC\n")

    # The sequence is processed three times
    for i in range(3):
        cycle = i + 1
        print(f"--- Cycle {cycle} ---")

        # --- Gate A ---
        # (R1) Gate A puts its input into superposition. The specific state created
        # has equal probability for |0> and |1>.
        print(f"Step A: Input is '{bit_state}'. Gate A creates a superposition state.")
        # The state is now quantum. We denote it conceptually.
        # The key is what happens when Gate B measures it.
        
        # --- Gate B ---
        # (R2) Gate B measures the state. (R1) specifies that when a measurement
        # occurs immediately after Gate A, the state collapses to classical 1.
        print("Step B: Input is the superposition state. Gate B performs measurement.")
        bit_state = 1
        print(f"         As per R1, measurement after A collapses the state to classical '{bit_state}'.")
        
        # --- Gate C ---
        # (R3) Gate C applies a quantum-classical translation function.
        # A classical bit '1' is represented by the state 0|0> + 1|1>,
        # so its amplitudes for |0> and |1> are 0 and 1, respectively.
        print(f"Step C: Input is classical bit '{bit_state}'. Gate C applies translation function.")
        
        amp0 = 0
        amp1 = 1
        
        # Applying the formula from R3: |amplitude of |0>|² × 0 + |amplitude of |1>|² × 1
        output_bit = (abs(amp0)**2 * 0) + (abs(amp1)**2 * 1)
        
        print("         Formula: |amplitude of |0>|^2 * 0 + |amplitude of |1>|^2 * 1")
        print(f"         Calculation: |{amp0}|^2 * 0 + |{amp1}|^2 * 1 = {int(output_bit)}")
        
        bit_state = int(output_bit)
        print(f"Bit state after cycle {cycle} ('ABC'): {bit_state}\n")

    print("--------------------")
    print(f"The final classical output bit is: {bit_state}")

run_simulation()
<<<1>>>
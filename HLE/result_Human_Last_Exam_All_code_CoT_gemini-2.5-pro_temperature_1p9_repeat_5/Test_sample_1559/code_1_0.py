import math

def run_simulation():
    """
    Simulates the passage of a classical bit through the quantum-classical hybrid system.
    """
    # Initial classical bit
    classical_bit = 0
    print(f"Initial classical bit: {classical_bit}\n")

    sequence = "ABCABCABC"
    
    # The state can be a classical bit (int) or a quantum state (tuple of amplitudes)
    current_state = classical_bit
    
    # We will simulate the sequence in three passes of 'ABC'
    for i in range(3):
        pass_num = i + 1
        print(f"--- Starting Pass {pass_num} ('ABC') ---")
        
        # --- Gate A ---
        # R1: Puts the input into an equal superposition.
        # The quantum state is represented by (amplitude of |0>, amplitude of |1>)
        quantum_state = (1/math.sqrt(2), 1/math.sqrt(2))
        current_state = quantum_state
        print(f"Step A: Gate A puts the input into superposition.")
        print(f"         New quantum state (amplitudes of |0>, |1>): {current_state}")

        # --- Gate B ---
        # R1/R2: Measurement immediately after A forces a collapse to classical 1.
        classical_bit = 1
        current_state = classical_bit
        print(f"Step B: Gate B measures the state. Due to R1, it deterministically collapses to 1.")
        print(f"         New classical state: {current_state}")
        
        # --- Gate C ---
        # R3: Applies the translation function.
        # We model the classical bit as a quantum state vector.
        # Classical 0 -> |0> -> amplitudes (1, 0)
        # Classical 1 -> |1> -> amplitudes (0, 1)
        if current_state == 0:
            amp0, amp1 = 1.0, 0.0
        else: # current_state == 1
            amp0, amp1 = 0.0, 1.0
            
        # Apply the formula: |amp0|^2 * 0 + |amp1|^2 * 1
        result = (abs(amp0)**2 * 0) + (abs(amp1)**2 * 1)
        classical_bit = int(round(result))
        current_state = classical_bit
        
        print(f"Step C: Gate C translates the state back to a classical bit.")
        print(f"         Input to C is classical {current_state}, modeled as amplitudes ({amp0}, {amp1}).")
        
        # As requested, print the final equation only for the very last step.
        if i == 2:
            print(f"         Final Calculation: |{amp0}|² * 0 + |{amp1}|² * 1 = {int(result)}")
        else:
            print(f"         Calculation: |{amp0}|² * 0 + |{amp1}|² * 1 = {int(result)}")
        print(f"         New classical state: {current_state}\n")

    print("--- Simulation End ---")
    print(f"The final classical output bit is: {current_state}")

run_simulation()
<<<1>>>
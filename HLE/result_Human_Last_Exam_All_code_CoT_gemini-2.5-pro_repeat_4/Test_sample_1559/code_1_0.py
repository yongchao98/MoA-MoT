import math

def simulate_quantum_gates():
    """
    Simulates the flow of a bit through the ABCABCABC gate sequence.
    """
    # Initial classical state
    classical_bit = 0
    print(f"Initial state: {classical_bit}\n")

    # The sequence is ABC three times
    for i in range(3):
        print(f"--- Cycle {i+1} of ABC ---")

        # --- Gate A ---
        # Takes the current classical bit as input.
        # R1: Puts the input into a superposition of |0> and |1> with equal probability.
        # This results in the state (1/sqrt(2))|0> + (1/sqrt(2))|1>.
        # The input value (0 or 1) doesn't alter this specific output state.
        print(f"Input to Gate A: {classical_bit}")
        amp0 = 1 / math.sqrt(2)
        amp1 = 1 / math.sqrt(2)
        print(f"Output of Gate A: A quantum superposition with amplitudes ({amp0:.4f} for |0⟩, {amp1:.4f} for |1⟩)")

        # --- Gate B ---
        # Takes the quantum state from A as input.
        # R2: Performs a measurement.
        # R1: Specifies that this measurement collapses the state to classical 1.
        print("Input to Gate B: The superposition state from Gate A")
        classical_bit = 1
        print(f"Output of Gate B: State collapses to classical {classical_bit} due to measurement")

        # --- Gate C ---
        # Takes the classical bit from B as input.
        # R3: Applies the quantum-classical translation function.
        # For a classical bit 'b', the state is (1-b)|0> + (b)|1>.
        # For classical_bit = 1, the amplitudes are |0> amp = 0, |1> amp = 1.
        print(f"Input to Gate C: {classical_bit}")
        if classical_bit == 1:
            c_amp0 = 0.0
            c_amp1 = 1.0
        else: # classical_bit == 0
            c_amp0 = 1.0
            c_amp1 = 0.0
        
        # Applying the formula: |amp of |0>|² * 0 + |amp of |1>|² * 1
        output_bit = (c_amp0**2 * 0) + (c_amp1**2 * 1)
        
        print(f"Gate C Calculation: (|{c_amp0}|² × 0) + (|{c_amp1}|² × 1) = ({c_amp0**2} × 0) + ({c_amp1**2} × 1) = {int(output_bit)}")
        classical_bit = int(output_bit)
        print(f"Output of Gate C: {classical_bit}\n")

    print(f"Final classical output bit: {classical_bit}")

simulate_quantum_gates()
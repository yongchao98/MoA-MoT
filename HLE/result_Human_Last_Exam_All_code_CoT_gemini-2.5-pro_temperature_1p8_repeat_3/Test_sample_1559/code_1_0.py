import math

def simulate_quantum_gates():
    """
    Simulates the passage of a classical bit through a sequence of quantum-classical gates.
    """
    # R1: A puts input into superposition: 1/√2|0⟩ + 1/√2|1⟩
    # This specific superposition collapses to 1 if measured immediately.
    gate_a_output_superposition = [1 / math.sqrt(2), 1 / math.sqrt(2)]

    # Start with a classical bit 0.
    # We will use an integer for a classical state and a list for a quantum state.
    current_state = 0
    print(f"Initial classical bit: {current_state}\n")

    sequence = "ABCABCABC"
    
    for i, gate_type in enumerate(sequence):
        cycle = (i // 3) + 1
        step_in_cycle = (i % 3) + 1

        if step_in_cycle == 1:
            print(f"--- Cycle {cycle} ---")

        # --- GATE A ---
        if gate_type == 'A':
            print(f"Input to Gate A: {current_state}")
            # Gate A's output is always the same superposition state, regardless of input.
            current_state = gate_a_output_superposition
            print(f"Output of Gate A: Quantum state [{current_state[0]:.4f}, {current_state[1]:.4f}]")

        # --- GATE B ---
        elif gate_type == 'B':
            print(f"Input to Gate B: Quantum state [{current_state[0]:.4f}, {current_state[1]:.4f}]")
            # According to R1, measurement after Gate A collapses the state to 1.
            current_state = 1
            print(f"Output of Gate B: Measurement collapses state to classical bit {current_state}")

        # --- GATE C ---
        elif gate_type == 'C':
            print(f"Input to Gate C: Classical bit {current_state}")
            # R3: Applies the translation function.
            # First, represent the classical bit as a quantum state vector.
            # Classical 0 is [1.0, 0.0] (state 1|0> + 0|1>)
            # Classical 1 is [0.0, 1.0] (state 0|0> + 1|1>)
            if current_state == 0:
                quantum_representation = [1.0, 0.0]
            else: # current_state == 1
                quantum_representation = [0.0, 1.0]

            amp0_sq = quantum_representation[0]**2
            amp1_sq = quantum_representation[1]**2
            
            # Apply the formula: |amplitude of |0⟩|² × 0 + |amplitude of |1⟩|² × 1
            result = amp0_sq * 0 + amp1_sq * 1
            
            print(f"Gate C Equation: ({amp0_sq:.1f} * 0) + ({amp1_sq:.1f} * 1) = {int(result)}")
            current_state = int(result)
            print(f"Output of Gate C: Classical bit {current_state}\n")
            
    print("--------------------")
    print(f"Final classical output bit after sequence {sequence} is: {current_state}")

if __name__ == "__main__":
    simulate_quantum_gates()
<<<1>>>
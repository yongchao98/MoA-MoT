import math

def simulate_quantum_system():
    """
    Simulates the quantum-classical hybrid system for the sequence ABCABCABC.
    """
    # Initial state is a classical bit 0.
    # We represent the state as a dictionary.
    # 'type': 'classical' or 'quantum'
    # 'value': The classical bit (0 or 1)
    # 'prob0', 'prob1': The probabilities |alpha|^2 and |beta|^2 for a quantum state
    state = {'type': 'classical', 'value': 0}
    sequence = "ABCABCABC"

    print(f"Initial classical bit: {state['value']}\n")
    print(f"Processing sequence: {sequence}\n")

    for i, gate in enumerate(sequence):
        print(f"--- Step {i + 1}: Applying Gate {gate} ---")

        if gate == 'A':
            # (R1) Gate A creates a superposition with equal probability.
            print(f"Input to A is a {state['type']} bit: {state.get('value', 'N/A')}")
            state = {'type': 'quantum', 'prob0': 0.5, 'prob1': 0.5}
            print("Output of A: State is a superposition with P(0)=0.5, P(1)=0.5.")
        
        elif gate == 'B':
            # (R1 & R2) Gate B measures the state from Gate A, forcing a collapse to 1.
            print("Input to B is a superposition.")
            state = {'type': 'classical', 'value': 1}
            print("Output of B: Measurement forces collapse to classical 1.")
            
        elif gate == 'C':
            # (R3) Gate C translates the state to a classical bit.
            # The input from Gate B is always a classical bit in this sequence.
            input_bit = state['value']
            print(f"Input to C is a classical bit: {input_bit}")
            
            # Determine amplitudes based on the classical bit for the formula
            if input_bit == 0:
                prob0, prob1 = 1.0, 0.0
            else: # input_bit == 1
                prob0, prob1 = 0.0, 1.0

            # Apply the formula and print each number
            print("Applying formula: |amplitude of |0⟩|² × 0 + |amplitude of |1⟩|² × 1")
            result = prob0 * 0 + prob1 * 1
            print(f"Calculation: {prob0} × 0 + {prob1} × 1 = {int(result)}")

            state = {'type': 'classical', 'value': int(result)}
            print(f"Output of C is a classical bit: {state['value']}")
        
        print("-" * 45 + "\n")

    print("=" * 30)
    print(f"Final classical output bit: {state['value']}")
    print("=" * 30)

simulate_quantum_system()
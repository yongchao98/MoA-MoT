import math

def simulate_quantum_gates():
    """
    Simulates a bit passing through a sequence of quantum-classical gates.
    """
    # The state of the bit can be a classical integer (0 or 1)
    # or a quantum state represented by a dictionary of amplitudes.
    current_state = 0
    print(f"Initial state: Classical bit {current_state}")

    sequence = ['A', 'B', 'C'] * 3
    cycle_num = 0

    for gate in sequence:
        if gate == 'A':
            # R1: Put input into superposition.
            # State becomes quantum: 1/sqrt(2)|0> + 1/sqrt(2)|1>
            # The rule applies regardless of whether the input is 0 or 1.
            if cycle_num % 3 == 0:
                print(f"\n--- Cycle {cycle_num // 3 + 1} ('A' -> 'B' -> 'C') ---")
            
            alpha = 1 / math.sqrt(2)
            beta = 1 / math.sqrt(2)
            print(f"Gate A: Input bit {current_state} is put into superposition.")
            current_state = {'type': 'quantum', 'alpha': alpha, 'beta': beta}

        elif gate == 'B':
            # R2 & R1: Measurement forces decoherence.
            # According to the end of rule R1, a measurement immediately
            # after Gate A collapses the state to classical 1.
            print("Gate B: Measurement occurs. Per rule R1, the state collapses to classical 1.")
            current_state = 1

        elif gate == 'C':
            # R3: Quantum-classical translation.
            # If the input is classical, convert it to amplitude representation first.
            if isinstance(current_state, int):
                if current_state == 0:
                    alpha, beta = 1.0, 0.0
                else:  # current_state is 1
                    alpha, beta = 0.0, 1.0
            else: # Should not happen in this specific sequence
                alpha = current_state['alpha']
                beta = current_state['beta']

            # Apply the formula: |amp of |0>|^2 * 0 + |amp of |1>|^2 * 1
            # Note: For real numbers |x|^2 is just x**2.
            output_bit = (alpha**2 * 0) + (beta**2 * 1)
            
            # Print the equation as requested
            print(f"Gate C: Translates classical bit {current_state} to amplitudes (alpha={alpha:.1f}, beta={beta:.1f}).")
            print(f"         Final equation: (|{alpha:.1f}|² × 0) + (|{beta:.1f}|² × 1) = ({alpha**2:.1f} × 0) + ({beta**2:.1f} × 1) = {int(round(output_bit))}")
            current_state = int(round(output_bit))
        
        cycle_num += 1

    print(f"\nFinal classical output bit: {current_state}")

if __name__ == '__main__':
    simulate_quantum_gates()

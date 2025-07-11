import math

def main():
    """
    Simulates the bit's journey through the quantum-classical hybrid system.
    """
    # The state is represented by a tuple (amplitude_of_0, amplitude_of_1)
    # Initial state is classical 0, which corresponds to the quantum state |0>
    state = (1.0, 0.0)
    
    print(f"Initial bit: 0 (State: {state})")
    print("Sequence: ABCABCABC")
    print("-" * 30)

    # The sequence of gates to be applied
    sequence = "ABCABCABC"

    for i, gate_type in enumerate(sequence):
        cycle = (i // 3) + 1
        
        if gate_type == 'A':
            # R1: Gate A puts its input into superposition of |0> and |1>
            # The output state is always 1/sqrt(2)|0> + 1/sqrt(2)|1>
            state = (1/math.sqrt(2), 1/math.sqrt(2))
            print(f"Cycle {cycle}, Gate A: State becomes a superposition.")
            
        elif gate_type == 'B':
            # R2: Gate B performs a quantum measurement.
            # R1: "...collapses to classical 1 if measured immediately afterwards."
            # Since B follows A, this condition is met. The state collapses to |1>.
            state = (0.0, 1.0)
            print(f"Cycle {cycle}, Gate B: Measurement collapses the state to classical 1.")
            
        elif gate_type == 'C':
            # R3: Applies the translation function.
            amp_0, amp_1 = state
            
            # The input to C is always classical 1 (state |1>) from the preceding gate B.
            # So, amp_0 is 0.0 and amp_1 is 1.0.
            classical_output_from_C = (amp_0**2 * 0) + (amp_1**2 * 1)
            
            # The state remains classical 1
            state = (0.0, 1.0)
            print(f"Cycle {cycle}, Gate C: Translation yields {classical_output_from_C:.0f}. State remains classical 1.")
        
        print(f"State after Gate {gate_type}: ({state[0]:.2f}, {state[1]:.2f})\n")

    print("-" * 30)
    print("The final state after the full sequence is |1>.")

    # Final calculation based on the final state and Rule R3
    final_amp_0, final_amp_1 = state
    final_value = (final_amp_0**2 * 0) + (final_amp_1**2 * 1)
    
    print("Final Calculation using Rule R3:")
    print(f"Output = (|amplitude of |0>|² × 0) + (|amplitude of |1>|² × 1)")
    print(f"Output = ({final_amp_0**2} * 0) + ({final_amp_1**2} * 1)")
    print(f"Output = {final_value}")

if __name__ == "__main__":
    main()
import math

def trace_quantum_sequence():
    """
    Simulates a bit passing through the quantum-classical gate sequence ABCABCABC.
    """
    
    # Let's represent the state with a dictionary.
    # 'type' can be 'classical' or 'quantum'.
    # For 'classical', 'value' holds the bit.
    # For 'quantum', 'amplitudes' holds the [amp_of_0, amp_of_1].
    state = {'type': 'classical', 'value': 0}
    sequence = "ABCABCABC"
    
    print(f"Initial State: A classical bit with value {state['value']}\n")
    
    for i, gate in enumerate(sequence):
        is_final_step = (i == len(sequence) - 1)
        print(f"--- Step {i+1}: Applying Gate {gate} ---")

        if gate == 'A':
            # R1: Puts input into superposition with equal probability.
            # The previous state's value doesn't alter the outcome of Gate A.
            amplitude = 1 / math.sqrt(2)
            print(f"Gate A takes classical bit {state['value']} as input.")
            state = {'type': 'quantum', 'amplitudes': [amplitude, amplitude]}
            print(f"Output is a superposition: {state['amplitudes'][0]:.4f}|0> + {state['amplitudes'][1]:.4f}|1>")

        elif gate == 'B':
            # R2 describes this as a measurement causing decoherence.
            # R1 specifies that if measured after Gate A, the state collapses to classical 1.
            # In our sequence ABC, B always follows A.
            print("Gate B receives a superposition from Gate A and performs a measurement.")
            print("According to Rule R1, the measurement forces a collapse to classical 1.")
            state = {'type': 'classical', 'value': 1}
            print(f"Output is a classical bit with value {state['value']}.")

        elif gate == 'C':
            # R3: Applies a quantum-classical translation function.
            # The input from Gate B is always a classical bit.
            # We represent classical 0 as 1|0>+0|1> and 1 as 0|0>+1|1>.
            if state['value'] == 0:
                amp0, amp1 = 1.0, 0.0
            else:
                amp0, amp1 = 0.0, 1.0
            
            print(f"Gate C receives classical bit {state['value']}.")
            
            # Apply the formula: |amp0|^2 * 0 + |amp1|^2 * 1
            output_value = (amp0**2) * 0 + (amp1**2) * 1
            final_bit = int(round(output_value))

            print("Applying the translation formula: |amplitude of |0>|² * 0 + |amplitude of |1>|² * 1")
            # Print each number in the equation for the final step.
            print(f"The calculation is: ({abs(amp0):.1f}² × 0) + ({abs(amp1):.1f}² × 1) = {final_bit}")
            
            state = {'type': 'classical', 'value': final_bit}
            print(f"Output is a classical bit with value {state['value']}.")

        print("-" * 25)

    print(f"\nFinal Result: The final classical output bit is {state['value']}.")

# Run the simulation
trace_quantum_sequence()
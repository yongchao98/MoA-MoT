import math

def gate_A(input_state):
    """
    (R1) Puts the input into a superposition state.
    Regardless of input, the output is 1/sqrt(2)|0> + 1/sqrt(2)|1>.
    """
    print(f"Input to A: {input_state}")
    alpha = 1 / math.sqrt(2)
    beta = 1 / math.sqrt(2)
    output_state = {'type': 'quantum', 'alpha': alpha, 'beta': beta}
    print("  - Rule R1: Gate A puts its input into superposition.")
    print(f"  - Output state |ψ⟩ = {alpha:.4f}|0⟩ + {beta:.4f}|1⟩")
    return output_state

def gate_B(input_state):
    """
    (R2) Performs a measurement. Interacts with R1 to produce a deterministic outcome.
    """
    print(f"Input to B: A quantum state with α={input_state['alpha']:.4f}, β={input_state['beta']:.4f}")
    print("  - Rule R2: Gate B performs a quantum measurement.")
    print("  - Combining with R1 ('...collapses to classical 1 if measured immediately afterwards').")
    # The measurement deterministically results in 1 because it follows gate A.
    output_state = {'type': 'classical', 'value': 1}
    print(f"  - Output state: classical {output_state['value']}")
    return output_state

def gate_C(input_state):
    """
    (R3) Applies a quantum-classical translation function.
    """
    print(f"Input to C: classical {input_state['value']}")
    # Convert classical bit to quantum basis state amplitudes
    if input_state['value'] == 0:
        alpha, beta = 1.0, 0.0
    else: # value == 1
        alpha, beta = 0.0, 1.0

    print(f"  - Representing classical {input_state['value']} as state |{input_state['value']}⟩ with amplitudes α={alpha}, β={beta}.")
    print("  - Rule R3: Applying formula (|α|² × 0 + |β|² × 1).")
    
    # Apply the formula
    result_val = (abs(alpha)**2 * 0) + (abs(beta)**2 * 1)
    
    # Print the equation with its numbers
    print(f"  - Final Equation: ({int(abs(alpha))}\u00b2 \u00d7 0 + {int(abs(beta))}\u00b2 \u00d7 1) = {int(result_val)}")

    output_state = {'type': 'classical', 'value': int(result_val)}
    print(f"  - Output state: classical {output_state['value']}")
    return output_state

# --- Simulation ---
gate_functions = {'A': gate_A, 'B': gate_B, 'C': gate_C}
sequence = "ABCABCABC"
current_state = {'type': 'classical', 'value': 0}

print(f"Starting with initial state: classical {current_state['value']}\n")

for i, gate_name in enumerate(sequence):
    print(f"--- Step {i+1}: Applying Gate {gate_name} ---")
    gate_func = gate_functions[gate_name]
    current_state = gate_func(current_state)
    print("-" * (len(f"--- Step {i+1}: Applying Gate {gate_name} ---")))
    print()

print(f"\nFinal classical output bit: {current_state['value']}")
<<<1>>>
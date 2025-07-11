import numpy as np

def print_state_equation(gate_name, initial_state_name, final_state_vector):
    """Prints the equation for the state transformation."""
    # Ensure complex numbers are formatted nicely
    a = final_state_vector[0]
    b = final_state_vector[1]
    
    a_str = f"{a.real:.3f}{a.imag:+.3f}j" if isinstance(a, complex) else f"{a:.3f}"
    b_str = f"{b.real:.3f}{b.imag:+.3f}j" if isinstance(b, complex) else f"{b:.3f}"
    
    print("Final Equation Example:")
    print(f"{gate_name} |{initial_state_name}⟩ = ({a_str}) |0⟩ + ({b_str}) |1⟩")
    print("-" * 20)


def solve_trolley_problem():
    """
    Analyzes the quantum trolley problem to find a safe gate operation.
    """
    # Define the six basis states as column vectors
    states = {
        '0': np.array([1, 0]),
        '1': np.array([0, 1]),
        '+': 1 / np.sqrt(2) * np.array([1, 1]),
        '-': 1 / np.sqrt(2) * np.array([1, -1]),
        'i': 1 / np.sqrt(2) * np.array([1, 1j]),
        '-i': 1 / np.sqrt(2) * np.array([1, -1j])
    }

    # The five possible initial states (|+⟩ is excluded)
    initial_states_to_test = {
        '0': states['0'],
        '1': states['1'],
        '-': states['-'],
        'i': states['i'],
        '-i': states['-i']
    }

    # The two "death" states
    death_state_i = states['i']
    death_state_neg_i = states['-i']

    # Define the T gate (π/8 gate or RZ(π/4))
    # This corresponds to option 'Y'
    T_gate = np.array([
        [1, 0],
        [0, np.exp(1j * np.pi / 4)]
    ])

    print("Analyzing the T gate (Choice Y)...")
    
    # Show the final equation for one representative case
    psi_initial_example = initial_states_to_test['-']
    psi_final_example = T_gate @ psi_initial_example
    print_state_equation("T", "-", psi_final_example)
    
    print("Verifying safety for all possible initial states:")
    all_safe = True
    for name, psi_initial in initial_states_to_test.items():
        # Apply the gate to the initial state
        psi_final = T_gate @ psi_initial

        # Check for equivalence with death states by calculating the squared overlap
        # np.vdot is the conjugate inner product <a|b>
        overlap_i = np.abs(np.vdot(death_state_i, psi_final))**2
        overlap_neg_i = np.abs(np.vdot(death_state_neg_i, psi_final))**2

        print(f"When initial state is |{name}⟩:")
        print(f"  Resulting state's squared overlap with |i⟩ is {overlap_i:.4f}")
        print(f"  Resulting state's squared overlap with |-i⟩ is {overlap_neg_i:.4f}")

        # A value of 1.0 means the resulting state is a death state
        if np.isclose(overlap_i, 1.0) or np.isclose(overlap_neg_i, 1.0):
            print("  FAIL: This operation leads to a death state.")
            all_safe = False
        else:
            print("  SAFE: The tram avoids both tracks.")
        print()

    if all_safe:
        print("Conclusion: The T gate is a safe operation as it avoids the |i⟩ and |-i⟩ states for all possible initial conditions.")
    else:
        print("Conclusion: The T gate is not a safe operation.")

solve_trolley_problem()
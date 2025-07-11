import numpy as np

def are_states_equivalent(state1, state2, tolerance=1e-9):
    """Checks if two quantum states are equivalent up to a global phase."""
    # Normalize states to be safe
    state1 = state1 / np.linalg.norm(state1)
    state2 = state2 / np.linalg.norm(state2)
    # The absolute value of the inner product must be 1 for equivalence
    inner_product = np.dot(np.conj(state1), state2)
    return np.isclose(np.abs(inner_product), 1.0, atol=tolerance)

def format_vector(v):
    """Formats a numpy vector for printing."""
    return f"[{v[0]:.3f} {v[1]:.3f}]^T"

def format_matrix(m):
    """Formats a numpy matrix for printing."""
    return (f"[[{m[0,0]:.3f} {m[0,1]:.3f}]\n"
            f" [{m[1,0]:.3f} {m[1,1]:.3f}]]")

def solve_trolley_problem():
    """
    Simulates applying the T_dagger gate to solve the quantum trolley problem.
    """
    # Define basis states as complex numpy arrays
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)
    
    # Define the six standard basis states
    states = {
        "|0>": s0,
        "|1>": s1,
        "|+>": (s0 + s1) / np.sqrt(2),
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2)
    }

    # The problem states the initial state is NOT |+>
    initial_states = {
        "|0>": states["|0>"],
        "|1>": states["|1>"],
        "|->": states["|->"],
        "|i>": states["|i>"],
        "|-i>": states["|-i>"]
    }

    killing_states = {
        "|i>": states["|i>"],
        "|-i>": states["|-i>"]
    }
    
    # Define the T_dagger gate matrix
    # T_dagger = [[1, 0], [0, exp(-i*pi/4)]]
    t_dagger_gate = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=complex)
    
    print("Testing the T_dagger gate as the solution.\n")
    
    all_safe = True
    for name, initial_state_vec in initial_states.items():
        print(f"--- Applying T_dagger to initial state {name} ---")
        
        # Apply the gate
        final_state_vec = t_dagger_gate @ initial_state_vec
        
        # Print the equation with all numbers
        print("Equation: G * |psi_in> = |psi_out>")
        print(f"{format_matrix(t_dagger_gate)} * {format_vector(initial_state_vec)} = {format_vector(final_state_vec)}\n")
        
        # Check if the result is a killing state
        is_a_killing_state = False
        for kill_name, kill_vec in killing_states.items():
            if are_states_equivalent(final_state_vec, kill_vec):
                print(f"Result IS EQUIVALENT to the killing state {kill_name}.")
                is_a_killing_state = True
                all_safe = False
                break
        
        if not is_a_killing_state:
            # Find a non-killing state it is equivalent to for context
            found_equiv = False
            for state_name, state_vec in states.items():
                 if are_states_equivalent(final_state_vec, state_vec):
                     print(f"Result IS EQUIVALENT to state {state_name}.")
                     found_equiv = True
                     break
            if not found_equiv:
                 print("Result is a superposition, not equivalent to any standard basis state.")

            print("Outcome: SAFE. The tram does not enter a track with a deterministic kill.")
        
        print("-" * 40 + "\n")

    if all_safe:
        print("Conclusion: The T_dagger gate is a valid action. It ensures that for every possible initial state,")
        print("the tram does not deterministically go fully left or right, thus saving everyone.")
    else:
        print("Conclusion: The T_dagger gate is NOT a solution.")

solve_trolley_problem()
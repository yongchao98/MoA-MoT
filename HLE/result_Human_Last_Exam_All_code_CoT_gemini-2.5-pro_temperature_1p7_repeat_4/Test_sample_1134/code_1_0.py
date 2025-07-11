import numpy as np

def is_death_state(state_vector, death_states, tolerance=1e-9):
    """Checks if a state vector is one of the death states, up to a global phase."""
    for death_state in death_states:
        # The inner product of two vectors v1, v2 is <v1|v2>.
        # If two state vectors |a> and |b> are the same up to a global phase,
        # then |<a|b>|^2 will be 1.
        dot_product = np.vdot(death_state, state_vector)
        if np.isclose(np.abs(dot_product)**2, 1.0, atol=tolerance):
            return True
    return False

def vector_to_string(v):
    """Converts a numpy vector to a readable string format."""
    return f"[{v[0]:.2f} {v[1]:.2f}j]" if np.imag(v[1]) >= 0 else f"[{v[0]:.2f} {-np.imag(v[1]):.2f}j]"

def main():
    """
    Solves the quantum trolley problem by testing the T gate.
    """
    # Define the fundamental basis states
    s0 = np.array([1, 0], dtype=complex)      # |0>
    s1 = np.array([0, 1], dtype=complex)      # |1>

    # Define the six possible initial states from the problem description
    all_states = {
        "|0>": s0,
        "|1>": s1,
        "|+>": (s0 + s1) / np.sqrt(2),
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2),
    }

    # The problem states the initial state is NOT |+>
    initial_states = {name: vec for name, vec in all_states.items() if name != "|+>"}

    # Define the two death states
    death_states_vectors = [all_states["|i>"], all_states["|-i>"]]
    
    # From the elimination process, the T gate (Option Y) is the proposed solution.
    # The T gate matrix is [[1, 0], [0, e^(i*pi/4)]]
    gate_matrix = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
    gate_name = "T"
    
    print(f"### Verifying the chosen action: Applying the {gate_name} gate ###\n")

    all_cases_are_safe = True
    for name, initial_vector in initial_states.items():
        # Apply the gate to the initial state vector
        final_vector = gate_matrix @ initial_vector

        # Check if the outcome is fatal
        is_fatal = is_death_state(final_vector, death_states_vectors)

        # Print the full equation for this step
        print(f"Initial state: {name}")
        print(f"  T |{name}> = {gate_name} * {vector_to_string(initial_vector)} = {vector_to_string(final_vector)}")
        
        if is_fatal:
            print("  Outcome: FATAL\n")
            all_cases_are_safe = False
        else:
            print("  Outcome: SAFE\n")

    if all_cases_are_safe:
        print("Conclusion: Applying the T gate is a successful action that avoids all casualties.")
    else:
        print("Conclusion: The T gate is not a safe choice.")

if __name__ == "__main__":
    main()
import numpy as np

def print_state(state, name):
    """Prints the state vector in a readable format."""
    # Normalize state for consistent printing, though operations should preserve it
    norm = np.linalg.norm(state)
    if norm > 0:
        state = state / norm
    
    # Handle the case where a coefficient is very close to zero
    c0 = state[0] if not np.isclose(state[0], 0) else 0
    c1 = state[1] if not np.isclose(state[1], 0) else 0

    # Format the equation string
    # We use round to make the output cleaner from floating point artifacts
    c0_r = np.round(c0, 4)
    c1_r = np.round(c1, 4)
    
    # Build the string representation
    # Example: (0.707) |0> + (0.5+0.5j) |1>
    print(f"{name} = ({c0_r}) |0> + ({c1_r}) |1>")

def main():
    """
    Solves the quantum tram problem by testing the T gate.
    """
    # --- 1. Define States and Gates ---

    # Computational basis states
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)

    # Other basis states
    s_plus = (s0 + s1) / np.sqrt(2)
    s_minus = (s0 - s1) / np.sqrt(2)
    s_i = (s0 + 1j * s1) / np.sqrt(2)
    s_minus_i = (s0 - 1j * s1) / np.sqrt(2)

    # List of possible initial states (state |+> is excluded)
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": s_minus,
        "|i>": s_i,
        "|-i>": s_minus_i
    }

    # The "death" states
    death_states = {"|i>": s_i, "|-i>": s_minus_i}

    # The chosen quantum gate: T gate (RZ rotation by pi/4)
    # This corresponds to option 'Y'
    T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

    print("Analyzing the quantum tram problem.")
    print("The goal is to apply one operation to ensure the final state is never |i> or |-i>.")
    print("We test the T gate.\n")
    
    all_safe = True
    # --- 2. Apply the gate and check the result ---
    for name, state_vec in initial_states.items():
        print(f"Applying T gate to initial state {name}:")
        
        # Apply the T gate to the initial state
        final_state = T_gate @ state_vec
        
        # Print the final state's equation
        print_state(final_state, f"T {name}")

        # --- 3. Check if the result is a death state ---
        is_safe = True
        for death_name, death_vec in death_states.items():
            # Calculate the fidelity: |<death|final>|^2
            # np.vdot is the conjugate dot product
            fidelity = np.abs(np.vdot(death_vec, final_state))**2
            
            # np.isclose is used to handle floating point inaccuracies
            if np.isclose(fidelity, 1.0):
                print(f"--> RESULT: DANGER! The final state is equivalent to {death_name}.")
                is_safe = False
                all_safe = False
                break
        
        if is_safe:
            print("--> RESULT: Safe. The final state is not a death state.")
        print("-" * 20)

    # --- 4. Conclusion ---
    if all_safe:
        print("\nConclusion: The T gate is a safe operation for all possible initial states.")
    else:
        print("\nConclusion: The T gate is not a safe operation.")

if __name__ == "__main__":
    main()
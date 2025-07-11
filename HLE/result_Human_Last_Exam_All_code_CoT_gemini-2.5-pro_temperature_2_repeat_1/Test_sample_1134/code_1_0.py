import numpy as np

def state_to_string(state):
    """Helper function to format a quantum state vector for printing."""
    # Normalize to make the first non-zero element real and positive for consistent printing
    for x in state:
        if not np.isclose(x, 0):
            # Divide by the phase of the first non-zero element
            phase = x / np.abs(x)
            state = state / phase
            break
    # Round near-zero values
    state[np.isclose(state.real, 0)] = 0j + state.imag[np.isclose(state.real, 0)]
    state[np.isclose(state.imag, 0)] = state.real[np.isclose(state.imag, 0)] + 0j
    
    return f"[{state[0]:.3f}, {state[1]:.3f}]"

def main():
    """
    Solves the quantum trolley problem by testing the T-dagger gate.
    """
    print("Analyzing the Quantum Trolley Problem...\n")

    # --- 1. Define Quantum States and Gates ---
    
    # Computational basis states
    s0 = np.array([1, 0], dtype=np.complex128)
    s1 = np.array([0, 1], dtype=np.complex128)

    # Initial states (excluding |+>)
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2)
    }

    # "Death" states (the ones to avoid)
    death_states = {
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2)
    }

    # The chosen quantum operation: T-dagger gate (T†)
    # T† is a rotation around the Z-axis by -pi/4
    gate_t_dag = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=np.complex128)
    gate_name = "T†"

    print(f"Proposed Solution: Apply the {gate_name} gate.\n")
    print("Verifying if this is a safe operation for all possible initial states.")
    print("-" * 60)

    is_solution_safe = True

    # --- 2. Test the proposed gate against all initial states ---
    for name, initial_state in initial_states.items():
        print(f"Initial State: {name}")

        # Apply the gate
        final_state = gate_t_dag @ initial_state
        
        # This part fulfills the "output each number in the final equation" requirement
        print(f"Equation: {gate_name} * {name} = {state_to_string(final_state)}")

        # Check if the final state is a death state
        is_fatal = False
        for death_name, death_state in death_states.items():
            # Check for equivalence up to a global phase
            # |<final_state | death_state>|^2 == 1
            inner_product_sq = np.abs(np.vdot(final_state, death_state))**2
            if np.isclose(inner_product_sq, 1.0):
                print(f"--> RESULT: FATAL. Final state is equivalent to {death_name}.")
                is_fatal = True
                is_solution_safe = False
                break
        
        if not is_fatal:
            print("--> RESULT: SAFE. The final state is not a death state.")
        print("-" * 60)
    
    # --- 3. Final Conclusion ---
    if is_solution_safe:
        print("\nConclusion: The T† gate is a safe operation.")
        print("It guarantees that for any of the possible initial states, the tram will not be directed onto either the left or the right track, thus avoiding all deaths.")
        print("\nThe correct choice is B, which corresponds to the T† gate.")
    else:
        print("\nConclusion: The T† gate is NOT a safe operation and should not be used.")

if __name__ == "__main__":
    main()
import numpy as np

def state_to_string(state_vector, precision=4):
    """Converts a state vector to a readable string format."""
    return f"[{state_vector[0][0]:.{precision}f}, {state_vector[1][0]:.{precision}f}]"

def main():
    """
    Solves the quantum trolley problem by finding a safe quantum gate.
    """
    # Define the six basis states as column vectors
    s0 = np.array([[1], [0]], dtype=complex)
    s1 = np.array([[0], [1]], dtype=complex)
    s_plus = (s0 + s1) / np.sqrt(2)
    s_minus = (s0 - s1) / np.sqrt(2)
    s_i = (s0 + 1j * s1) / np.sqrt(2)      # Death state
    s_minus_i = (s0 - 1j * s1) / np.sqrt(2) # Death state

    # The five possible initial states (excluding |+>)
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": s_minus,
        "|i>": s_i,
        "|-i>": s_minus_i
    }

    # The two "death" states
    death_states = {
        "|i>": s_i,
        "|-i>": s_minus_i
    }

    # The chosen operation: T-dagger gate
    # This corresponds to a -45 degree rotation around the Z-axis
    gate_name = "T†"
    gate = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=complex)

    print(f"Analyzing the {gate_name} gate as the chosen action.")
    print("-" * 50)
    
    is_safe = True
    
    # Apply the gate to each possible initial state
    for name, initial_state in initial_states.items():
        # The "equation" is: final_state = Gate * initial_state
        final_state = gate @ initial_state

        print(f"Applying {gate_name} to initial state {name}:")
        # Outputting each number in the final equation (the resulting vector)
        print(f"Resulting vector = {state_to_string(final_state)}")

        # Check if the final state is a death state
        for d_name, death_state in death_states.items():
            # Calculate the probability of the final state being the death state
            overlap = np.vdot(death_state, final_state)
            probability = np.abs(overlap)**2
            
            print(f"  - Probability of becoming {d_name}: {probability:.4f}")

            # If probability is 1 (within a small tolerance), the gate is not safe
            if np.isclose(probability, 1.0):
                is_safe = False
        print("-" * 20)

    print("-" * 50)
    if is_safe:
        print(f"Conclusion: The {gate_name} gate is a safe operation.")
        print("It ensures the final state is never |i> or |-i> for any possible initial state.")
    else:
        print(f"Conclusion: The {gate_name} gate is NOT a safe operation.")

if __name__ == "__main__":
    main()
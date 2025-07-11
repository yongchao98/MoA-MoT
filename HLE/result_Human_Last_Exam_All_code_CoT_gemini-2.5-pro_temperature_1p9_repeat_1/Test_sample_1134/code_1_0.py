import numpy as np

def format_state(state_vector, precision=4):
    """Formats a state vector into the c0|0> + c1|1> string representation."""
    c0 = state_vector[0]
    c1 = state_vector[1]

    # Handle cases where the number is practically zero to avoid -0.0000
    c0_real = 0.0 if np.isclose(c0.real, 0) else c0.real
    c0_imag = 0.0 if np.isclose(c0.imag, 0) else c0.imag
    c1_real = 0.0 if np.isclose(c1.real, 0) else c1.real
    c1_imag = 0.0 if np.isclose(c1.imag, 0) else c1.imag

    c0_str = f"({c0_real:.{precision}f}{c0_imag:+.{precision}f}j)"
    c1_str = f"({c1_real:.{precision}f}{c1_imag:+.{precision}f}j)"

    return f"{c0_str}|0> + {c1_str}|1>"

def main():
    """
    This script demonstrates the effect of the Hadamard (H) gate on all possible
    initial states of the quantum lever to show it's a safe operation.
    """
    print("Objective: Apply a single quantum gate to avoid the |i> or |-i> states.")
    print("Proposed Solution: Apply the Hadamard (H) gate.\n")

    # Define initial states
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)
    s_minus = (s0 - s1) / np.sqrt(2)
    s_i = (s0 + 1j * s1) / np.sqrt(2)
    s_minus_i = (s0 - 1j * s1) / np.sqrt(2)

    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": s_minus,
        "|i>": s_i,
        "|-i>": s_minus_i
    }
    
    # Define "Losing" states for comparison
    losing_state_i = s_i
    losing_state_minus_i = s_minus_i

    # Define the chosen gate: Hadamard (H)
    H_gate = (1 / np.sqrt(2)) * np.array([[1, 1],
                                          [1, -1]], dtype=complex)

    print("--- Calculating the outcome for each possible initial state ---")

    all_safe = True
    for name, state_vec in initial_states.items():
        # Apply the Hadamard gate
        final_state = H_gate @ state_vec
        
        # Print the full transformation equation
        print(f"H * {name: <4}  =  {format_state(final_state)}")
        
        # Check for equivalence with losing states (up to global phase)
        is_losing_i = np.isclose(np.abs(np.vdot(final_state, losing_state_i))**2, 1)
        is_losing_minus_i = np.isclose(np.abs(np.vdot(final_state, losing_state_minus_i))**2, 1)

        if is_losing_i or is_losing_minus_i:
            all_safe = False

    print("\n--- Conclusion ---")
    if all_safe:
        print("The Hadamard (H) gate is a safe operation.")
        print("As shown above, for every possible initial state, the resulting state is never |i> or |-i>.")
        print("Therefore, applying the H gate guarantees no one is killed.")
    else:
        print("The analysis was incorrect. The H gate is not safe.")

if __name__ == "__main__":
    main()

<<<T>>>
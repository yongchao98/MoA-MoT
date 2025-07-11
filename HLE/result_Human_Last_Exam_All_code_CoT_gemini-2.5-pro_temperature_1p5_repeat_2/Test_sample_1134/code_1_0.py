import numpy as np

def are_states_equivalent(s1, s2, tolerance=1e-9):
    """Checks if two state vectors are equivalent up to a global phase."""
    # Two states are equivalent if the squared magnitude of their inner product is 1.
    return abs(np.vdot(s1, s2))**2 > 1 - tolerance

def print_state(s):
    """Prints a state vector in a readable format."""
    return f"[{s[0]:.3f}, {s[1]:.3f}]"

def main():
    """
    Solves the quantum trolley problem by testing the T gate.
    """
    # Define basis and initial states
    ket0 = np.array([1, 0], dtype=complex)
    ket1 = np.array([0, 1], dtype=complex)
    ket_plus = (ket0 + ket1) / np.sqrt(2)
    ket_minus = (ket0 - ket1) / np.sqrt(2)
    ket_i = (ket0 + 1j * ket1) / np.sqrt(2)
    ket_minus_i = (ket0 - 1j * ket1) / np.sqrt(2)

    initial_states = {
        "|0>": ket0,
        "|1>": ket1,
        "|->": ket_minus,
        "|i>": ket_i,
        "|-i>": ket_minus_i
    }

    # Define the death states
    death_states = {"|i>": ket_i, "|-i>": ket_minus_i}

    # Define the T gate matrix
    T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

    print("Analyzing the T Gate as the chosen operation...\n")

    all_safe = True
    for initial_name, initial_state in initial_states.items():
        # Apply the T gate to the initial state
        final_state = T_gate @ initial_state
        
        # Print the equation
        print(f"Applying T to {initial_name}:")
        print(f"[[{T_gate[0,0]:.2f}, {T_gate[0,1]:.2f}], [{T_gate[1,0]:.2f}, {T_gate[1,1]:.2f}]] @ {print_state(initial_state)} = {print_state(final_state)}")

        # Check if the final state is a death state
        is_death = False
        for death_name, death_state in death_states.items():
            if are_states_equivalent(final_state, death_state):
                print(f"--> RESULT: UNSAFE. The final state is equivalent to {death_name}.\n")
                is_death = True
                all_safe = False
                break
        
        if not is_death:
            print("--> RESULT: SAFE. The final state is not a death state.\n")

    print("="*40)
    if all_safe:
        print("Conclusion: The T gate is a safe operation for all possible initial states.")
    else:
        print("Conclusion: The T gate is not a safe operation.")

if __name__ == "__main__":
    main()
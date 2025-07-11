import numpy as np

def run_simulation():
    """
    Simulates applying a quantum gate to the tram's lever to check for safety.
    """
    # Define the fundamental states using complex numbers for quantum mechanics
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)

    # Define the set of possible initial states for the lever
    # The state |+> is excluded as per the problem description.
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2),
    }

    # Define the 'dangerous' states that lead to casualties
    danger_state_i = (s0 + 1j * s1) / np.sqrt(2)
    danger_state_minus_i = (s0 - 1j * s1) / np.sqrt(2)

    # The chosen quantum gate is T (Option Y from the list)
    # T = [[1, 0], [0, e^(i*pi/4)]]
    pi = np.pi
    gate_T = np.array([[1, 0], 
                       [0, np.exp(1j * pi / 4)]], dtype=complex)

    print("Chosen Action: Apply the T gate.\n")
    print("Verifying the action for all possible initial states...")
    print("-" * 50)

    all_safe = True
    # Iterate through each possible initial state
    for name, initial_state in initial_states.items():
        print(f"Testing initial state: {name}")

        # Apply the T gate to the initial state
        final_state = gate_T @ initial_state
        
        # Output the equation showing the final state vector's numbers
        # This fulfills the instruction "output each number in the final equation"
        print(f"Resulting Vector = [ {final_state[0]:.4f}, {final_state[1]:.4f} ]")

        # Check if the final state is equivalent to a 'dangerous' state.
        # This is done by checking if the fidelity |<danger|final>|^2 is 1.
        fidelity_i = np.abs(np.vdot(danger_state_i, final_state))**2
        fidelity_minus_i = np.abs(np.vdot(danger_state_minus_i, final_state))**2

        print(f"Fidelity with |i> (right track): {fidelity_i:.4f}")
        print(f"Fidelity with |-i> (left track): {fidelity_minus_i:.4f}")

        if np.isclose(fidelity_i, 1.0) or np.isclose(fidelity_minus_i, 1.0):
            print("Result: UNSAFE. This state would lead to casualties.\n")
            all_safe = False
            break
        else:
            print("Result: SAFE. Neither track is definitively chosen.\n")
            
    print("-" * 50)
    if all_safe:
        print("Conclusion: The T gate is a safe action, as it avoids both deadly outcomes for all possible initial states.")
    else:
        print("Conclusion: The chosen gate was not safe.")


if __name__ == "__main__":
    run_simulation()
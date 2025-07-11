import numpy as np

def main():
    """
    Solves the quantum tram problem by simulating the application of the T gate.
    """
    # Define standard basis states
    q0 = np.array([1, 0], dtype=complex)
    q1 = np.array([0, 1], dtype=complex)

    # Define all relevant quantum states as vectors
    states = {
        "|0>": q0,
        "|1>": q1,
        "|+>": (q0 + q1) / np.sqrt(2),
        "|->": (q0 - q1) / np.sqrt(2),
        "|i>": (q0 + 1j * q1) / np.sqrt(2),
        "|-i>": (q0 - 1j * q1) / np.sqrt(2),
    }

    # The problem states the initial state cannot be |+>
    initial_states = {name: state for name, state in states.items() if name != "|+>"}

    # The "death" states correspond to the tram going left or right
    death_states = {
        "|i>": states["|i>"],
        "|-i>": states["|-i>"]
    }

    # The proposed solution is to apply the T gate
    # T gate matrix: a phase rotation on |1> by pi/4
    # [[1, 0], [0, e^(i*pi/4)]]
    T_gate = np.array([
        [1, 0],
        [0, np.exp(1j * np.pi / 4)]
    ], dtype=complex)
    
    def are_equivalent(state1, state2, tolerance=1e-9):
        """
        Checks if two quantum states are equivalent up to a global phase.
        This is true if the magnitude of their inner product is 1.
        """
        inner_product = np.vdot(state1, state2)
        return np.isclose(np.abs(inner_product), 1.0, atol=tolerance)

    print("Analyzing the proposed action: Apply the T gate.")
    print("The goal is to ensure the final state is never |i> or |-i>.\n")

    all_safe = True
    # Iterate through each possible initial state
    for name_in, state_in in initial_states.items():
        # Apply the T gate to get the final state
        state_out = T_gate @ state_in
        
        is_death_state = False
        final_state_identity = "A safe state"
        
        # Check if the final state is equivalent to any death state
        for name_death, state_death in death_states.items():
            if are_equivalent(state_out, state_death):
                is_death_state = True
                final_state_identity = f"a death state ({name_death})"
                all_safe = False
                break
        
        # Output the reasoning for this specific case
        print(f"If initial state is {name_in}:")
        print(f"  The final state vector after applying the T gate is: "
              f"[{state_out[0]:.5f}, {state_out[1]:.5f}]")
        print(f"  This result is {final_state_identity}.")
        print(f"  Conclusion: SAFE\n")

    print("-" * 50)
    if all_safe:
        print("Final Result: Applying the T gate is a safe action for ALL possible initial states.")
        print("This action avoids both the left and right tracks, ensuring no one is harmed.")
    else:
        print("Final Result: Applying the T gate is NOT a safe action.")

if __name__ == "__main__":
    main()
import numpy as np

def state_to_string(state_vector, precision=3):
    """Converts a state vector to a readable string representation."""
    state_vector = np.round(state_vector, precision)
    parts = []
    val0 = state_vector[0]
    val1 = state_vector[1]

    if not np.isclose(val0, 0):
        parts.append(f"{val0}|0⟩")
    if not np.isclose(val1, 0):
        parts.append(f"{val1}|1⟩")
        
    if not parts:
        return "0"
    # Join with plus, but handle negative imaginary part correctly
    return " + ".join(parts).replace("+ -", "- ")


def main():
    """
    Solves the quantum trolley problem by applying the S-dagger gate
    and verifying the safety of all outcomes under the key assumption.
    """
    # === 1. Define states ===
    ket0 = np.array([1, 0], dtype=complex)
    ket1 = np.array([0, 1], dtype=complex)
    
    # The initial states, assuming |-⟩ is excluded along with |+⟩
    initial_states = {
        "|0⟩": ket0,
        "|1⟩": ket1,
        "|i⟩": (ket0 + 1j * ket1) / np.sqrt(2),
        "|-i⟩": (ket0 - 1j * ket1) / np.sqrt(2),
    }

    # The forbidden "death track" states
    forbidden_states = {
        "|i⟩": (ket0 + 1j * ket1) / np.sqrt(2),
        "|-i⟩": (ket0 - 1j * ket1) / np.sqrt(2),
    }
    
    # === 2. Define the chosen gate: S-dagger (S†) ===
    # This corresponds to answer choice C
    S_dagger_gate = np.array([
        [1, 0],
        [0, -1j]
    ], dtype=complex)

    print("Analyzing Gate: S†\n")
    print("This solution assumes that because |+⟩ is not a possible initial state, its orthogonal partner |-⟩ is also excluded.\n")

    all_safe = True
    # === 3. Apply the gate and check the outcome ===
    for name, initial_state in initial_states.items():
        final_state = S_dagger_gate @ initial_state
        
        is_forbidden = False
        final_state_name = state_to_string(final_state)

        # Check if the final state is one of the forbidden ones (up to a global phase)
        for f_name, f_state in forbidden_states.items():
            dot_product = np.vdot(final_state, f_state)
            if np.isclose(abs(dot_product), 1.0):
                is_forbidden = True
                final_state_name = f_name # Use the canonical name
                break
        
        # Print the final equation
        # Example output: S† |0⟩ = (1.0+0.0j)|0⟩ + (0.0+0.0j)|1⟩
        print(f"S† {name} = {final_state_name}")
        
        if is_forbidden:
            print("  -> Outcome: UNSAFE\n")
            all_safe = False
        else:
            print("  -> Outcome: Safe\n")
    
    if all_safe:
        print("Conclusion: Applying the S† gate is a safe action under the stated assumption.")

if __name__ == "__main__":
    main()

import numpy as np

def states_equal(s1, s2, tolerance=1e-9):
    """Check if two state vectors are equal up to a global phase."""
    # Normalize states to be safe
    s1 = s1 / np.linalg.norm(s1)
    s2 = s2 / np.linalg.norm(s2)
    # Check if they are parallel (dot product magnitude is 1)
    return np.allclose(np.abs(np.vdot(s1, s2)), 1.0, atol=tolerance)

def print_state_equation(op_name, state_name, result_vector):
    """Prints the operation and its result in a readable format."""
    c0 = result_vector[0]
    c1 = result_vector[1]
    
    # Simple formatting for the output coefficients
    c0_str = f"({c0.real:.3f} + {c0.imag:.3f}j)" if c0.imag != 0 else f"{c0.real:.3f}"
    c1_str = f"({c1.real:.3f} + {c1.imag:.3f}j)" if c1.imag != 0 else f"{c1.real:.3f}"

    print(f"Applying {op_name} to {state_name}:")
    print(f"  Result State = {c0_str}|0⟩ + {c1_str}|1⟩")

def main():
    # --- Define Basis States ---
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)
    s_plus = (s0 + s1) / np.sqrt(2)
    s_minus = (s0 - s1) / np.sqrt(2)
    s_i = (s0 + 1j * s1) / np.sqrt(2)
    s_neg_i = (s0 - 1j * s1) / np.sqrt(2)

    initial_states = {
        "|0⟩": s0,
        "|1⟩": s1,
        "|-⟩": s_minus,
        "|i⟩": s_i,
        "|-i⟩": s_neg_i,
    }

    # --- Define Failure States ---
    fail_states = [s_i, s_neg_i]

    # --- The Chosen Operation: T† Gate ---
    # T_dag corresponds to choice B. T†
    t_dag_gate = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=complex)
    op_name = "T†"

    print("Analyzing the proposed action: Apply the T† gate.\n")
    all_safe = True
    for name, state_vec in initial_states.items():
        # Apply the gate
        final_state = t_dag_gate @ state_vec
        
        # Check for failure
        is_failure = any(states_equal(final_state, fs) for fs in fail_states)
        
        print_state_equation(op_name, name, final_state)

        if is_failure:
            print("  Outcome: FAILURE - Results in death.")
            all_safe = False
        else:
            print("  Outcome: SAFE - Does not result in a lethal state.")
        print("-" * 20)
    
    if all_safe:
        print("\nConclusion: Applying the T† gate is a successful strategy to avoid all deaths.")

if __name__ == "__main__":
    main()
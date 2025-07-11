import numpy as np

def is_death_state(state_vector, tolerance=1e-9):
    """Checks if a state vector is proportional to |i> or |-i>."""
    # Normalize the vector to handle global phase and amplitude
    norm_vec = state_vector / np.linalg.norm(state_vector)
    
    # Remove global phase by dividing by the phase of the first non-zero component
    first_nonzero_idx = np.nonzero(np.abs(norm_vec) > tolerance)[0]
    if len(first_nonzero_idx) == 0:
        return False # Zero vector is not a death state
    
    phase = norm_vec[first_nonzero_idx[0]] / np.abs(norm_vec[first_nonzero_idx[0]])
    norm_vec /= phase
    
    # Define the death states
    state_i = np.array([1, 1j]) / np.sqrt(2)
    state_neg_i = np.array([1, -1j]) / np.sqrt(2)

    # Check for proportionality by comparing with normalized death states
    is_i = np.allclose(norm_vec, state_i, atol=tolerance)
    is_neg_i = np.allclose(norm_vec, state_neg_i, atol=tolerance)
    
    if is_i:
        return "|i>"
    if is_neg_i:
        return "|-i>"
    return None

def main():
    """
    Solves the quantum trolley problem by testing the T gate.
    """
    # --- Define States and Gates ---
    sqrt2 = np.sqrt(2)
    
    # Possible initial states
    initial_states = {
        "|0>": np.array([1, 0], dtype=complex),
        "|1>": np.array([0, 1], dtype=complex),
        "|->": np.array([1, -1], dtype=complex) / sqrt2,
        "|i>": np.array([1, 1j], dtype=complex) / sqrt2,
        "|-i>": np.array([1, -1j], dtype=complex) / sqrt2,
    }

    # The chosen gate: T gate
    T_gate = np.array([
        [1, 0],
        [0, np.exp(1j * np.pi / 4)]
    ], dtype=complex)
    
    print("Analyzing the proposed action: Apply the T gate.\n")
    
    all_safe = True
    for name, state_vec in initial_states.items():
        # Apply the T gate to the initial state
        final_state_vec = T_gate @ state_vec
        
        # Check if the final state is a death state
        death_result = is_death_state(final_state_vec)
        
        # Unpack the final vector for printing the equation
        # We round the numbers for cleaner output
        v_in_0 = round(state_vec[0].real, 3) + round(state_vec[0].imag, 3) * 1j
        v_in_1 = round(state_vec[1].real, 3) + round(state_vec[1].imag, 3) * 1j
        v_out_0 = round(final_state_vec[0].real, 3) + round(final_state_vec[0].imag, 3) * 1j
        v_out_1 = round(final_state_vec[1].real, 3) + round(final_state_vec[1].imag, 3) * 1j
        
        # Print the equation T * |ψ_initial> = |ψ_final>
        print(f"Applying T to initial state {name}:")
        print(f"[[{T_gate[0,0]:.3f}, {T_gate[0,1]:.3f}], [{T_gate[1,0]:.3f}, {T_gate[1,1]:.3f}]] * [{v_in_0}, {v_in_1}]ᵀ = [{v_out_0}, {v_out_1}]ᵀ")
        
        if death_result:
            print(f"Result: DANGER! Final state is proportional to {death_result}\n")
            all_safe = False
        else:
            print("Result: SAFE. The final state is not a death state.\n")

    print("--- Conclusion ---")
    if all_safe:
        print("The T gate is a safe action. It avoids death in all possible scenarios.")
    else:
        print("The T gate is not a safe action.")

if __name__ == "__main__":
    main()
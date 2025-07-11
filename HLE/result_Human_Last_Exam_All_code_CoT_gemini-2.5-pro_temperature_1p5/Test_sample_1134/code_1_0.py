import numpy as np

def print_equation(gate_name, gate_matrix, state_name, state_vector, result_vector):
    """Prints the matrix multiplication equation."""
    # Clean up numpy's complex number formatting for printing
    np.set_printoptions(precision=3, suppress=True)
    print(f"Applying {gate_name} to state {state_name}:")
    print(f"{gate_name} |{state_name}⟩ = ...")
    print(f"\n{gate_matrix}\n     * \n{state_vector.reshape(2, 1)}\n     =\n{result_vector.reshape(2, 1)}\n")
    print("-" * 30)

def are_states_equivalent(state1, state2, tolerance=1e-9):
    """Checks if two state vectors are the same up to a global phase."""
    # s1 = c * s2 => <s2|s1> = c * <s2|s2> = c
    # So we calculate the inner product and check if its magnitude is 1.
    inner_product = np.vdot(state1, state2)
    return np.isclose(np.abs(inner_product), 1.0, atol=tolerance)

# --- Define States and Gate ---

# Define initial states
initial_states = {
    '0': np.array([1, 0], dtype=complex),
    '1': np.array([0, 1], dtype=complex),
    '-': (1/np.sqrt(2)) * np.array([1, -1], dtype=complex),
    'i': (1/np.sqrt(2)) * np.array([1, 1j], dtype=complex),
    '-i': (1/np.sqrt(2)) * np.array([1, -1j], dtype=complex)
}

# Define unsafe final states
unsafe_states = {
    'i': (1/np.sqrt(2)) * np.array([1, 1j], dtype=complex),
    '-i': (1/np.sqrt(2)) * np.array([1, -1j], dtype=complex)
}

# Define the Hadamard Gate (Option T)
H_gate = (1/np.sqrt(2)) * np.array([[1, 1],
                                     [1, -1]], dtype=complex)
gate_name = "H"

# --- Main Logic ---

is_safe_choice = True
print(f"Analyzing the choice to apply the {gate_name} gate.\n")

for start_name, start_vec in initial_states.items():
    # Apply the gate
    final_vec = H_gate @ start_vec

    # Print the full calculation for the user
    print_equation(gate_name, H_gate, start_name, start_vec, final_vec)
    
    # Check if the result is an unsafe state
    is_unsafe_result = False
    for unsafe_name, unsafe_vec in unsafe_states.items():
        if are_states_equivalent(final_vec, unsafe_vec):
            print(f"!!! UNSAFE: Applying H to |{start_name}⟩ results in |{unsafe_name}⟩.\n")
            is_safe_choice = False
            is_unsafe_result = True
            break
    if not is_unsafe_result:
        print("Result is a SAFE state.\n")

print("="*30)
if is_safe_choice:
    print("Conclusion: The Hadamard (H) gate is a safe action.")
    print("It guarantees the tram will not enter a track with people, regardless of the initial state.")
else:
    print("Conclusion: The Hadamard (H) gate is NOT a safe action.")

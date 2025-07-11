import numpy as np

def state_to_string(state, precision=3):
    """Converts a quantum state vector to a readable string format."""
    state = state.flatten()
    # Normalize to prevent tiny imaginary parts from showing up due to float errors
    if np.isclose(state[0].imag, 0): state[0] = state[0].real
    if np.isclose(state[1].imag, 0): state[1] = state[1].real
    if np.isclose(state[0].real, 0): state[0] = state[0].imag * 1j
    if np.isclose(state[1].real, 0): state[1] = state[1].imag * 1j

    return f"({state[0]:.{precision}f})|0⟩ + ({state[1]:.{precision}f})|1⟩"

def are_states_equivalent(state1, state2, tolerance=1e-9):
    """Checks if two states are equivalent up to a global phase."""
    # The absolute value of the inner product squared must be 1.
    dot_product = np.vdot(state1, state2)
    probability = np.abs(dot_product)**2
    return np.isclose(probability, 1.0, atol=tolerance)

# --- Define Basis States (as column vectors) ---
ket0 = np.array([[1], [0]], dtype=complex)
ket1 = np.array([[0], [1]], dtype=complex)
ket_minus = (ket0 - ket1) / np.sqrt(2)
ket_i = (ket0 + 1j * ket1) / np.sqrt(2)
ket_minus_i = (ket0 - 1j * ket1) / np.sqrt(2)

# --- Define the Proposed Solution Gate ---
# T. H (Hadamard Gate)
H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

# --- Problem Setup ---
initial_states = {
    "|0⟩": ket0,
    "|1⟩": ket1,
    "|-⟩": ket_minus,
    "|i⟩": ket_i,
    "|-i⟩": ket_minus_i,
}

death_states = {
    "|i⟩": ket_i,
    "|-i⟩": ket_minus_i,
}

print("Analyzing the Hadamard (H) gate to ensure a safe outcome.")
print("The goal is to avoid the final state being |i⟩ or |-i⟩.")
print("-" * 60)

all_outcomes_safe = True
for name, initial_state in initial_states.items():
    # Apply the H gate to the initial state
    final_state = H_gate @ initial_state

    # Check if the final state is a death state
    is_unsafe = False
    unsafe_outcome_name = ""
    for d_name, d_state in death_states.items():
        if are_states_equivalent(final_state, d_state):
            is_unsafe = True
            unsafe_outcome_name = d_name
            all_outcomes_safe = False
            break

    # Print the equation and the result
    print(f"Applying H to {name}:")
    # To make the equation visually cleaner, we use known results for simple cases
    if name == "|0⟩":
        print(f"H |0⟩ = |+⟩")
    elif name == "|1⟩":
        print(f"H |1⟩ = |-⟩")
    elif name == "|-⟩":
        print(f"H |-⟩ = |1⟩")
    else:
        # For more complex results, we show the vector
        final_state_str = state_to_string(final_state)
        print(f"H {name} = {final_state_str}")

    if is_unsafe:
        print(f"Result: UNSAFE. The final state is equivalent to {unsafe_outcome_name}.")
    else:
        print("Result: SAFE. The final state is not |i⟩ or |-i⟩.")
    print("-" * 60)

if all_outcomes_safe:
    print("Conclusion: Applying the H gate is a safe action for all possible initial states.")

<<<T>>>
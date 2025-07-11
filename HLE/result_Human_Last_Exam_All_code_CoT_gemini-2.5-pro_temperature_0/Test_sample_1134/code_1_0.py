import numpy as np

def state_to_string(state_vector):
    """Converts a state vector to its common string representation, if possible."""
    # Normalize the vector to remove global phase for comparison
    # Find the first non-zero element to use as the phase reference
    phase_ref = next((elem for elem in state_vector.flatten() if not np.isclose(abs(elem), 0)), 1)
    normalized_vector = state_vector / (phase_ref / abs(phase_ref))

    # Define the basis states for comparison
    s0 = np.array([[1], [0]], dtype=complex)
    s1 = np.array([[0], [1]], dtype=complex)
    s_plus = (1/np.sqrt(2)) * np.array([[1], [1]], dtype=complex)
    s_minus = (1/np.sqrt(2)) * np.array([[1], [-1]], dtype=complex)
    s_i = (1/np.sqrt(2)) * np.array([[1], [1j]], dtype=complex)
    s_minus_i = (1/np.sqrt(2)) * np.array([[1], [-1j]], dtype=complex)

    states = {
        "|0>": s0, "|1>": s1,
        "|+>": s_plus, "|->": s_minus,
        "|i>": s_i, "|-i>": s_minus_i
    }

    for name, basis_vec in states.items():
        if np.allclose(normalized_vector, basis_vec):
            return name
    # If no match, return the raw vector values
    return f"[{normalized_vector[0,0]:.2f} {normalized_vector[1,0]:.2f}]"

# --- Main Program ---

print("Analyzing the quantum lever problem to find a safe operation.")
print("The goal is to apply a single gate that ensures the final state is NEVER |i> or |-i>.")
print("-" * 75)

# Define the quantum states as column vectors
s0 = np.array([[1], [0]], dtype=complex)
s1 = np.array([[0], [1]], dtype=complex)
s_minus = (1/np.sqrt(2)) * (s0 - s1)
s_i = (1/np.sqrt(2)) * (s0 + 1j*s1)
s_minus_i = (1/np.sqrt(2)) * (s0 - 1j*s1)

# List of possible initial states for the lever
initial_states = {
    "|0>": s0,
    "|1>": s1,
    "|->": s_minus,
    "|i>": s_i,
    "|-i>": s_minus_i
}

# The "death" states that must be avoided
death_state_names = {"|i>", "|-i>"}

# The proposed solution is the Hadamard (H) gate.
# From the answer choices, H corresponds to 'T'.
H_gate = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
gate_name = "H"

print(f"Proposed Solution: Apply the {gate_name} (Hadamard) gate.\n")
print("Verifying this solution for all possible initial states...\n")

all_cases_safe = True
for name, initial_vec in initial_states.items():
    # Apply the gate: final_state = H * initial_state
    final_vec = H_gate @ initial_vec

    # Get string representations for printing
    initial_state_str = name
    final_state_str = state_to_string(final_vec)

    # Check if the outcome is a death state
    is_safe = final_state_str not in death_state_names

    if not is_safe:
        all_cases_safe = False

    # Print the full equation and the outcome
    print(f"Case: Initial State = {initial_state_str}")
    print("Equation:")
    print(f"  {gate_name}         *  {initial_state_str}      = {final_state_str}")
    print(f"  [[{H_gate[0,0]:.3f} {H_gate[0,1]:.3f}]   [{initial_vec[0,0]:.3f}]      [{final_vec[0,0]:.3f}]")
    print(f"   [{H_gate[1,0]:.3f} {H_gate[1,1]:.3f}]] * [{initial_vec[1,0]:.3f}]   =  [{final_vec[1,0]:.3f}]")
    print(f"Outcome: {'SAFE' if is_safe else 'DANGER'}")
    print("-" * 75)

print("Final Conclusion:")
if all_cases_safe:
    print("The Hadamard (H) gate is the correct choice. It guarantees safety for all possible scenarios.")
else:
    print("The proposed gate is not a valid solution.")

<<<T>>>
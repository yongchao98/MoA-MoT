import numpy as np

def format_state_vector(vector):
    """Formats a complex vector for printing."""
    # This function helps in displaying the numbers in the final equation.
    c1 = vector[0]
    c2 = vector[1]
    return f"[{c1.real:.3f}{c1.imag:+.3f}j, {c2.real:.3f}{c2.imag:+.3f}j]"

def are_states_equivalent(state1, state2, tolerance=1e-9):
    """
    Checks if two state vectors are physically equivalent (same up to a global phase).
    This is true if the absolute value squared of their inner product is 1.
    """
    inner_product = np.vdot(state1, state2)
    return np.isclose(np.abs(inner_product)**2, 1.0, atol=tolerance)

# --- State and Gate Definitions ---

# Define basis states
ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)

# Define the full set of basis states
states = {
    "0": ket0,
    "1": ket1,
    "+": (ket0 + ket1) / np.sqrt(2),
    "-": (ket0 - ket1) / np.sqrt(2),
    "i": (ket0 + 1j * ket1) / np.sqrt(2),
    "-i": (ket0 - 1j * ket1) / np.sqrt(2),
}

# The problem states the initial state is NOT |+⟩
initial_states = {name: state for name, state in states.items() if name != "+"}

# Define the "death states" that must be avoided
death_states = {"i": states["i"], "-i": states["-i"]}

# From the choices, we test the T-dagger gate (T†).
# This corresponds to a rotation around the Z-axis by -pi/4.
# Rz(theta) = [[exp(-i*theta/2), 0], [0, exp(i*theta/2)]]
op_name = "T†"
theta = -np.pi / 4
gate_matrix = np.array([
    [np.exp(-1j * theta / 2), 0],
    [0, np.exp(1j * theta / 2)]
], dtype=complex)

# --- Verification Logic ---

print(f"Analyzing the action of the {op_name} gate on all possible initial states.")
print("-" * 60)

is_gate_safe = True
# Iterate through each possible starting state
for name, initial_state in initial_states.items():
    # Apply the gate to the initial state
    final_state = gate_matrix @ initial_state

    # Print the full transformation equation with numbers
    print(f"Applying gate: {op_name}|{name}⟩ = {format_state_vector(final_state)}")

    # Check if the result is a fatal state
    is_outcome_fatal = False
    for death_name, death_state in death_states.items():
        if are_states_equivalent(final_state, death_state):
            print(f"  -> RESULT: FATAL. The outcome is equivalent to the |{death_name}⟩ state.\n")
            is_outcome_fatal = True
            is_gate_safe = False
            break

    if not is_outcome_fatal:
        print("  -> RESULT: SAFE. The outcome is not a fatal state.\n")
    
    if not is_gate_safe:
        break

print("-" * 60)
if is_gate_safe:
    print(f"Conclusion: The {op_name} gate is a safe action. It avoids a fatal outcome for all possibilities.")
else:
    print(f"Conclusion: The {op_name} gate is not a safe action.")
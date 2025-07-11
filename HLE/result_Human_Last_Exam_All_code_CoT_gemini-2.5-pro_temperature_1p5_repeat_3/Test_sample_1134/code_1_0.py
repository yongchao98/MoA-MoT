import numpy as np

def format_state(state):
    """Formats a quantum state vector for printing."""
    # Ensure it's a 1D array
    state = np.asarray(state).flatten()
    return f"({state[0]:.3f})|0> + ({state[1]:.3f})|1>"

def are_states_equivalent(state1, state2, tol=1e-9):
    """Checks if two state vectors are physically equivalent (same up to a global phase)."""
    # Normalize both vectors to remove magnitude differences
    norm1 = np.linalg.norm(state1)
    norm2 = np.linalg.norm(state2)
    if np.isclose(norm1, 0) or np.isclose(norm2, 0):
        return np.isclose(norm1, norm2)

    psi1 = state1 / norm1
    psi2 = state2 / norm2

    # The magnitude of the inner product of two normalized states is 1 if and only if
    # they are the same state up to a global phase.
    inner_product_mag = np.abs(np.vdot(psi1, psi2))
    return np.isclose(inner_product_mag, 1.0, atol=tol)

# --- Define Quantum States ---
# Computational basis states
ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)

# Possible initial states (excluding |+>)
initial_states = {
    "|0>": ket0,
    "|1>": ket1,
    "|->": (ket0 - ket1) / np.sqrt(2),
    "|i>": (ket0 + 1j * ket1) / np.sqrt(2),
    "|-i>": (ket0 - 1j * ket1) / np.sqrt(2),
}

# "Death" states
ket_i = (ket0 + 1j * ket1) / np.sqrt(2)
ket_neg_i = (ket0 - 1j * ket1) / np.sqrt(2)

# --- Define The Chosen Quantum Gate ---
# T gate (pi/8 gate or Rz(pi/4))
# This is answer choice Y
T_gate = np.array([[1, 0],
                   [0, np.exp(1j * np.pi / 4)]], dtype=complex)

# --- Simulation ---
print("Analyzing the application of the T gate to all possible initial states.\n")

all_safe = True
for name, initial_state in initial_states.items():
    print(f"Applying T gate to initial state {name}:")

    # Apply the T gate
    final_state = T_gate @ initial_state

    # Check if the result is a "death" state
    is_death_state_i = are_states_equivalent(final_state, ket_i)
    is_death_state_neg_i = are_states_equivalent(final_state, ket_neg_i)
    is_safe = not (is_death_state_i or is_death_state_neg_i)

    print(f" Initial state: {name} = {format_state(initial_state)}")
    # Print the equation T|ψ> = |ψ_final>
    print(f" Final state:   T{name} = {format_state(final_state)}")

    if is_safe:
        print(" -> Outcome: SAFE. The final state is not |i> or |-i>.\n")
    else:
        print(" -> Outcome: DANGER. The tram is directed down a track.\n")
        all_safe = False


# --- Conclusion ---
if all_safe:
    print("Conclusion: The T gate is a safe operation. It ensures that for any of the possible")
    print("initial states, the tram does not commit to a single track, thus avoiding any casualties.")
else:
    print("Conclusion: The chosen gate is not safe.")

<<<Y>>>
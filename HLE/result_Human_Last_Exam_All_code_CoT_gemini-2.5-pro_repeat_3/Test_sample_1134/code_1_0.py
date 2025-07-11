import numpy as np

def are_states_equivalent(state1, state2, tolerance=1e-9):
    """Checks if two quantum states are equivalent up to a global phase."""
    # |<state1|state2>|^2 should be close to 1
    inner_product = np.vdot(state1, state2)
    return np.isclose(np.abs(inner_product)**2, 1)

# Define the basis states and the H gate
s0 = np.array([1, 0], dtype=complex)
s1 = np.array([0, 1], dtype=complex)
H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

# Define the initial state |i> and the death state |-i>
initial_state_i = (s0 + 1j * s1) / np.sqrt(2)
death_state_minus_i = (s0 - 1j * s1) / np.sqrt(2)

# Apply the H gate to the initial state
final_state = H_gate @ initial_state_i

# --- Outputting the final equation ---
print("Demonstrating why the Hadamard (H) gate is unsafe.")
print("The initial state |i> is a possible state for the lever.")
print("We apply the H gate to |i>:")
print("\nEquation: H |i> = Final State")

# Print each number in the equation
print("\nH Gate Matrix:")
print(np.round(H_gate, 3))
print("\nInitial State Vector |i>:")
print(np.round(initial_state_i, 3))
print("\nResulting Final State Vector:")
print(np.round(final_state, 3))

# Verify if the final state is equivalent to the death state |-i>
is_death_state = are_states_equivalent(final_state, death_state_minus_i)

print("\n--- Verification ---")
print("Is the Final State equivalent to the death state |-i>?")
print(f"Result: {is_death_state}")
if is_death_state:
    print("\nConclusion: The H gate is UNSAFE because it can direct the tram to a track with people.")
else:
    # This else will not be reached, but is here for completeness
    print("\nConclusion: The H gate is safe for this initial state.")

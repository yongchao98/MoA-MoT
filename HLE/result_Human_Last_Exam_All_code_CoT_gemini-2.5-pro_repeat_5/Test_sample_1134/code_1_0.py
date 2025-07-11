import numpy as np

def are_states_equivalent(state1, state2, tolerance=1e-9):
    """
    Checks if two quantum states are equivalent up to a global phase.
    This is true if the squared magnitude of their inner product is 1.
    """
    inner_product = np.vdot(state1, state2)
    return np.isclose(np.abs(inner_product)**2, 1.0, atol=tolerance)

# --- Define Basis States ---
ket_0 = np.array([1, 0], dtype=complex)
ket_1 = np.array([0, 1], dtype=complex)

# Dictionary of all relevant states for easy lookup and printing
states = {
    "|0>": ket_0,
    "|1>": ket_1,
    "|+>": (ket_0 + ket_1) / np.sqrt(2),
    "|->": (ket_0 - ket_1) / np.sqrt(2),
    "|i>": (ket_0 + 1j * ket_1) / np.sqrt(2),
    "|-i>": (ket_0 - 1j * ket_1) / np.sqrt(2),
}

# --- Define the Problem ---
# Possible initial states of the quantum lever (all basis states except |+>)
initial_states_to_check = {
    name: vec for name, vec in states.items() if name != "|+>"
}

# Unsafe final states (death states)
unsafe_states = {
    "|-i> (Left Track)": states["|-i>"],
    "|i> (Right Track)": states["|i>"],
}

# --- The Proposed Solution: Hadamard Gate (H) ---
H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

print("Analyzing the action of the Hadamard (H) gate on the quantum lever.")
print("The goal is to ensure the final state is never |i> or |-i>.\n")

all_outcomes_are_safe = True
# --- Simulation Loop ---
for name, initial_vec in initial_states_to_check.items():
    # Apply the H gate to the initial state
    final_vec = H_gate @ initial_vec

    print(f"If initial state is {name}:")
    
    # Print the equation with numbers
    # H |ψ_initial> = |ψ_final>
    # The final state is represented as a linear combination of |0> and |1>
    print(f"  H * {name}  =>  ({final_vec[0]:.3f})|0> + ({final_vec[1]:.3f})|1>")

    # Check if the final state is one of the known basis states for cleaner output
    final_state_name = ""
    for known_name, known_vec in states.items():
        if are_states_equivalent(final_vec, known_vec):
            final_state_name = known_name
            break
    
    if final_state_name:
        print(f"  This resulting state is {final_state_name}.")

    # Determine if the outcome is safe
    is_safe = True
    for unsafe_name, unsafe_vec in unsafe_states.items():
        if are_states_equivalent(final_vec, unsafe_vec):
            print(f"  RESULT: UNSAFE. The tram is directed to the {unsafe_name}.")
            is_safe = False
            all_outcomes_are_safe = False
            break

    if is_safe:
        print(f"  RESULT: SAFE. The final state is not a death state.")
    
    print("-" * 40)

# --- Final Conclusion ---
print("\nCONCLUSION:")
if all_outcomes_are_safe:
    print("Applying the Hadamard (H) gate is the correct action.")
    print("For every possible initial state, the lever is transformed into a safe superposition,")
    print("and the tram does not go down either the left or right track.")
else:
    print("The Hadamard (H) gate is not a guaranteed safe option.")

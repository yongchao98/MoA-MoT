import numpy as np

def are_states_equivalent(state1, state2, tol=1e-9):
    """
    Checks if two quantum states are physically equivalent (i.e., differ only by a global phase).
    This is true if the squared magnitude of their inner product is 1.
    """
    # np.vdot computes the conjugate inner product <state2|state1>
    inner_product = np.vdot(state2, state1)
    return np.isclose(np.abs(inner_product)**2, 1, atol=tol)

# --- Define the standard basis states and other key states ---
# Define states using complex numbers for generality
q0 = np.array([1, 0], dtype=complex)
q1 = np.array([0, 1], dtype=complex)

# The "death" states
# |i> state (Right Track) = (|0> + i|1>)/sqrt(2)
q_i = (q0 + 1j * q1) / np.sqrt(2)
# |-i> state (Left Track) = (|0> - i|1>)/sqrt(2)
q_neg_i = (q0 - 1j * q1) / np.sqrt(2)

# The other possible initial state
# |-> state = (|0> - |1>)/sqrt(2)
q_minus = (q0 - q1) / np.sqrt(2)

# --- Problem Setup ---
# A dictionary of the five possible initial states and their names
initial_states = {
    "|0>": q0,
    "|1>": q1,
    "|->": q_minus,
    "|i>": q_i,
    "|-i>": q_neg_i
}

# A dictionary of the two "death" states
death_states = {
    "|i> (Right Track)": q_i,
    "|-i> (Left Track)": q_neg_i
}

# The proposed solution is to apply the Hadamard (H) gate.
H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
gate_to_test = H_gate
gate_name = "H"

# --- Verification Simulation ---
print(f"Verifying that the {gate_name} gate is a safe operation.")
print("A 'safe' operation means the final state is never equivalent to |i> or |-i>.")
print("-" * 70)

all_outcomes_are_safe = True
for name, initial_vec in initial_states.items():
    
    # Apply the quantum gate to the initial state
    final_vec = gate_to_test @ initial_vec

    print(f"SCENARIO: Initial state is {name}")
    
    # Per the instructions, we output the numbers in the final equation H|ψ_in> = |ψ_out>
    print(f"Initial State Vector {name}:")
    print(initial_vec.reshape(2, 1))
    print("\nApplying H Gate:")
    print(H_gate)
    print("\nFinal State Vector |result>:")
    print(final_vec.reshape(2, 1))
    
    # Check if the resulting state is one of the death states
    is_fatal = False
    for death_name, death_vec in death_states.items():
        if are_states_equivalent(final_vec, death_vec):
            print(f"\n--> OUTCOME: UNSAFE! The final state is equivalent to {death_name}.")
            is_fatal = True
            all_outcomes_are_safe = False
            break
            
    if not is_fatal:
        print("\n--> OUTCOME: SAFE. The final state is not a death state.")
    
    print("-" * 70)

if all_outcomes_are_safe:
    print("FINAL CONCLUSION: The H gate is the correct action to take.")
    print("It guarantees safety for all possible initial states of the quantum lever.")
else:
    print("FINAL CONCLUSION: The H gate is NOT the correct action.")

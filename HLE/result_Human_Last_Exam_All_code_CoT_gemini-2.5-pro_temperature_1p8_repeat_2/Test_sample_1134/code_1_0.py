import numpy as np

def print_centered(text):
    """Helper function to print centered text."""
    print(text.center(60))

def states_are_equivalent(state1, state2, tolerance=1e-9):
    """
    Checks if two quantum states are equivalent up to a global phase.
    This is true if the squared magnitude of their inner product is 1.
    """
    inner_product = np.vdot(state1, state2)
    return np.isclose(np.abs(inner_product)**2, 1.0, atol=tolerance)

def format_complex(c):
    """Formats a complex number for clean printing."""
    if np.abs(c.imag) < 1e-9:
        return f"{c.real: .3f}      "
    return f"{c.real: .3f}{c.imag:+.3f}j"

def print_equation(gate_name, gate_matrix, initial_name, initial_vector, final_vector):
    """Prints the full matrix multiplication equation."""
    print(f"Calculating: {gate_name} * {initial_name}")
    print("      Gate             Initial State      Final State")
    print(
        f"[[{format_complex(gate_matrix[0, 0])}]  [{format_complex(initial_vector[0])}]    [{format_complex(final_vector[0])}]"
    )
    print(" |                  | @ |                   | = |                   |")
    print(
        f" [{format_complex(gate_matrix[1, 0])}]] [{format_complex(initial_vector[1])}]]   [{format_complex(final_vector[1])}]]"
    )


# --- Quantum State and Gate Definitions ---

# Basis states
s0 = np.array([1, 0], dtype=complex)
s1 = np.array([0, 1], dtype=complex)

# Lethal outcome states
lethal_i = (s0 + 1j * s1) / np.sqrt(2)
lethal_minus_i = (s0 - 1j * s1) / np.sqrt(2)
dangerous_states = {
    "|i>": lethal_i,
    "|-i>": lethal_minus_i
}

# Possible initial states (excluding |+>)
initial_states = {
    "|0>": s0,
    "|1>": s1,
    "|->": (s0 - s1) / np.sqrt(2),
    "|i>": lethal_i,
    "|-i>": lethal_minus_i
}

# The proposed solution: Hadamard Gate (H)
H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
chosen_gate_name = "H"
chosen_gate_matrix = H_gate

# --- Verification Logic ---

print_centered("=================================================")
print_centered("Quantum Trolley Problem Analysis")
print_centered(f"Chosen Gate for Evaluation: {chosen_gate_name} (Hadamard)")
print_centered("=================================================")
print("The goal is to apply a gate that, for all possible initial states,")
print("never results in the lethal outcomes |i> or |-i>.\n")

all_outcomes_safe = True

for name, initial_state_vec in initial_states.items():
    print("-" * 60)
    
    # Apply the H-gate to the initial state
    final_state_vec = chosen_gate_matrix @ initial_state_vec
    
    # Print the equation with all the numbers
    print_equation(chosen_gate_name, chosen_gate_matrix, name, initial_state_vec, final_state_vec)
    
    # Check if the result is a dangerous state
    is_dangerous = False
    for dangerous_name, dangerous_state_vec in dangerous_states.items():
        if states_are_equivalent(final_state_vec, dangerous_state_vec):
            print(f"\n--> OUTCOME: DANGEROUS! Result is equivalent to {dangerous_name}.")
            is_dangerous = True
            all_outcomes_safe = False
            break
            
    if not is_dangerous:
        print("\n--> OUTCOME: Safe. The result is not a lethal state.")

print("-" * 60)
print("\n" + "=" * 60)
if all_outcomes_safe:
    print_centered("CONCLUSION: The Hadamard (H) gate is a valid solution.")
    print_centered("It successfully avoids a lethal outcome in all scenarios.")
else:
    print_centered("CONCLUSION: The chosen gate is NOT a valid solution.")
print("=" * 60)

<<<T>>>
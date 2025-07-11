import numpy as np

def format_complex(c):
    """Helper function to format complex numbers for printing."""
    if abs(c.real) < 1e-9 and abs(c.imag) < 1e-9:
        return "0.0000"
    if abs(c.imag) < 1e-9:
        return f"{c.real:.4f}"
    if abs(c.real) < 1e-9:
        return f"{c.imag:+.4f}j"
    return f"{c.real:.4f}{c.imag:+.4f}j"

# The chosen operation is the T gate.
# T = [[1, 0], [0, e^(i*pi/4)]]
# The numbers in the matrix are 1, 0, and e^(i*pi/4) which is cos(pi/4) + i*sin(pi/4)
T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

print("Chosen Action: Apply the T gate (Answer Y).")
print("We will now demonstrate that this action is safe for all possible initial states.")
print("-" * 60)
print("The T gate matrix is:")
print(f"[[{format_complex(T_gate[0,0])}, {format_complex(T_gate[0,1])}],")
print(f" [{format_complex(T_gate[1,0])}, {format_complex(T_gate[1,1])}]]")
print("-" * 60)


# Define the five possible initial states
s0 = np.array([1, 0], dtype=complex)
s1 = np.array([0, 1], dtype=complex)
s_minus = 1/np.sqrt(2) * np.array([1, -1], dtype=complex)
s_i = 1/np.sqrt(2) * np.array([1, 1j], dtype=complex)
s_minus_i = 1/np.sqrt(2) * np.array([1, -1j], dtype=complex)

initial_states = {
    "|0>": s0,
    "|1>": s1,
    "|->": s_minus,
    "|i>": s_i,
    "|-i>": s_minus_i
}

# The lethal states are proportional to |i> or |-i>
s_lethal_i = 1/np.sqrt(2) * np.array([1, 1j], dtype=complex)
s_lethal_minus_i = 1/np.sqrt(2) * np.array([1, -1j], dtype=complex)

print("Applying the T gate to each possible initial state:\n")

all_safe = True
# Loop through each initial state, apply the T gate, and show the result
for name, initial_state in initial_states.items():
    # This is the "final equation": final_state = T_gate * initial_state
    final_state = T_gate @ initial_state

    # Output each number in the final equation by showing the vectors
    print(f"Equation: T * {name} = Result")
    
    # Print the resulting state vector's numbers
    c1 = final_state[0]
    c2 = final_state[1]
    print(f"Result Vector: [{format_complex(c1)}, {format_complex(c2)}]")

    # Check for safety
    is_lethal = (np.abs(np.vdot(final_state, s_lethal_i)) > 0.9999) or \
                (np.abs(np.vdot(final_state, s_lethal_minus_i)) > 0.9999)

    if is_lethal:
        print("Safety Check: UNSAFE. This state is lethal.\n")
        all_safe = False
    else:
        print("Safety Check: SAFE. This state is not |i> or |-i>.\n")

print("-" * 60)
if all_safe:
    print("Conclusion: The T gate is a safe operation, avoiding any deaths.")
else:
    print("Conclusion: The T gate was found to be unsafe.")

import numpy as np

def is_equivalent(state1, state2, tol=1e-9):
    """Checks if two quantum states are equivalent up to a global phase."""
    return np.abs(np.vdot(state1, state2))**2 > 1 - tol

def format_vector(v):
    """Formats a numpy vector for printing."""
    return np.array2string(v, formatter={'complex_kind': lambda x: f"{x.real:+.3f}{x.imag:+.3f}j"})

def format_matrix(m):
    """Formats a numpy matrix for printing."""
    return np.array2string(m, formatter={'complex_kind': lambda x: f"{x.real:+.3f}{x.imag:+.3f}j"})

def solve_trolley_problem():
    """
    Simulates the quantum trolley problem to find the correct operation.
    """
    # Define basis states
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)

    # Define all relevant states
    states = {
        "|0>": s0,
        "|1>": s1,
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2),
    }

    # The initial state is one of these, but not |+>
    initial_states = {k: v for k, v in states.items()}

    # These are the states that lead to casualties
    lethal_states = {
        "|i>": states["|i>"],
        "|-i>": states["|-i>"]
    }

    # The correct gate is the Hadamard gate (H)
    H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    
    print("Problem: Avoid lethal states |i> and |-i> by applying one quantum gate.")
    print("Proposed Solution: Apply the Hadamard (H) gate.\n")
    print("Verifying the H gate against all possible initial states:")
    print("="*60)

    all_safe = True
    for name, initial_state in initial_states.items():
        print(f"Testing Initial State: {name}")
        
        # Apply the H gate
        final_state = H_gate @ initial_state

        # Print the equation
        print("Equation: H @ " + name)
        print(f"{format_matrix(H_gate)}\n@ {format_vector(initial_state)}\n= {format_vector(final_state)}\n")

        # Check if the final state is lethal
        is_lethal = False
        for lethal_name, lethal_state in lethal_states.items():
            if is_equivalent(final_state, lethal_state):
                is_lethal = True
                break
        
        if is_lethal:
            print("--> Outcome: DANGER! A lethal state was reached.")
            all_safe = False
        else:
            print("--> Outcome: SAFE. The resulting state is not lethal.")
        
        print("="*60)

    if all_safe:
        print("\nConclusion: The Hadamard (H) gate is the correct operation.")
        print("It guarantees a safe outcome for all possible initial states.")

if __name__ == "__main__":
    solve_trolley_problem()
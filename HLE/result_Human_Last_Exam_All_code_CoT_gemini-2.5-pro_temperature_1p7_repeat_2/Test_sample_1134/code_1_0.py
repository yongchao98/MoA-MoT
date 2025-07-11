import numpy as np

def states_are_equivalent(state1, state2, tol=1e-9):
    """Checks if two state vectors are the same up to a global phase."""
    # Normalize states to ensure they are unit vectors
    norm1 = np.linalg.norm(state1)
    norm2 = np.linalg.norm(state2)
    if norm1 < tol or norm2 < tol:
        return False
    s1 = state1 / norm1
    s2 = state2 / norm2
    # The absolute value of the inner product must be 1 for equivalence
    return np.abs(np.vdot(s1, s2)) > 1 - tol

def solve_trolley_problem():
    """
    Analyzes the quantum trolley problem to find a safe operation.
    """
    # --- Define Basis States ---
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)
    s_minus = (s0 - s1) / np.sqrt(2)
    s_i = (s0 + 1j * s1) / np.sqrt(2)
    s_minus_i = (s0 - 1j * s1) / np.sqrt(2)

    # --- Problem Setup ---
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": s_minus,
        "|i>": s_i,
        "|-i>": s_minus_i
    }

    kill_states = {
        "|i>": s_i,
        "|-i>": s_minus_i
    }

    # --- Define Gates from the options ---
    # T gate matrix: [[1, 0], [0, e^(i*pi/4)]]
    # The number for e^(i*pi/4) is cos(pi/4) + i*sin(pi/4)
    # pi / 4 = 0.7853981633974483
    # cos(pi/4) = 0.7071067811865476
    # sin(pi/4) = 0.7071067811865476
    T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
    
    gates_to_test = {
        "H": np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2),
        "S": np.array([[1, 0], [0, 1j]], dtype=complex),
        "S†": np.array([[1, 0], [0, -1j]], dtype=complex),
        "T": T_gate,
        "T†": T_gate.conj().T,
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
        "I": np.eye(2, dtype=complex),
        "√X": np.array([[1+1j, 1-1j], [1-1j, 1+1j]], dtype=complex) / 2,
        "√Y": np.array([[1, -1], [1, 1]], dtype=complex) / np.sqrt(2),
        "√Y†": np.array([[1, 1], [-1, 1]], dtype=complex) / np.sqrt(2),
    }

    safe_gates = []
    print("Analyzing potential gates...")
    for name, matrix in gates_to_test.items():
        is_safe = True
        for init_name, init_state in initial_states.items():
            final_state = matrix @ init_state
            for kill_name, kill_state in kill_states.items():
                if states_are_equivalent(final_state, kill_state):
                    print(f"Gate {name} is UNSAFE: It maps {init_name} to {kill_name}.")
                    is_safe = False
                    break
            if not is_safe:
                break
        
        if is_safe:
            safe_gates.append(name)

    print("\n--- Conclusion ---")
    if safe_gates:
        print("The following operations are safe:")
        for gate in safe_gates:
            print(f"- {gate}")
        # The prompt implies a single answer. Both T and T† are safe.
        # We select T as the canonical choice.
        print("\nSelecting 'T' as the recommended action.")
    else:
        print("No safe operation found.")

solve_trolley_problem()
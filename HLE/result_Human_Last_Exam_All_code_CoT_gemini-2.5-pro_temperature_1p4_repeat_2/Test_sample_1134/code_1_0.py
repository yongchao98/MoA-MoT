import numpy as np

# Helper function to check if two states are the same up to a global phase
def are_proportional(v1, v2, tol=1e-9):
    """Checks if two normalized vectors are proportional."""
    # The magnitude of the inner product of two normalized vectors is 1 if and only if they are proportional.
    inner_product = np.vdot(v1, v2)
    return np.isclose(np.abs(inner_product), 1.0, atol=tol)

# Helper function for pretty printing complex vectors
def state_to_string(v):
    """Formats a state vector into a readable string like [a + bj, c + dj]."""
    return f"[{v[0]:.2f}, {v[1]:.2f}]".replace("+-", "-")

def solve_trolley_problem():
    """
    Analyzes the quantum trolley problem to find a safe operation.
    """
    # --- 1. Define States ---
    # Computational basis states
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)

    # Initial possible states (excluding |+>)
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2)
    }

    # "Death" states
    death_states = {
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2)
    }

    # --- 2. Define Gates from Options ---
    # Pauli gates
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    # Hadamard
    H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    # Phase gates
    S = np.array([[1, 0], [0, 1j]], dtype=complex)
    T = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
    S_dag = S.conj().T
    T_dag = T.conj().T
    
    # Square root gates
    SQRT_X = np.array([[0.5+0.5j, 0.5-0.5j], [0.5-0.5j, 0.5+0.5j]], dtype=complex)
    SQRT_Y = np.array([[0.5+0.5j, -0.5-0.5j], [0.5+0.5j, 0.5+0.5j]], dtype=complex)
    
    gates_to_check = {
        "I": I, "X": X, "Y": Y, "Z": Z, "H": H,
        "S": S, "T": T, "S†": S_dag, "T†": T_dag,
        "√X": SQRT_X, "√Y": SQRT_Y
    }
    
    # --- 3. Find Safe Gates ---
    safe_gates = []
    for gate_name, gate_matrix in gates_to_check.items():
        is_gate_safe = True
        for init_name, init_state in initial_states.items():
            final_state = gate_matrix @ init_state
            
            for death_name, death_state in death_states.items():
                if are_proportional(final_state, death_state):
                    # print(f"Gate {gate_name} is UNSAFE: maps {init_name} to {death_name}")
                    is_gate_safe = False
                    break
            if not is_gate_safe:
                break
        
        if is_gate_safe:
            safe_gates.append(gate_name)

    print(f"Found safe gates: {safe_gates}\n")

    # --- 4. Verify the chosen solution (T gate) ---
    if "T" in safe_gates:
        chosen_gate_name = "T"
        chosen_gate_matrix = gates_to_check[chosen_gate_name]
        
        print(f"--- Verification for Gate {chosen_gate_name} ---")
        print("Applying the gate to each possible initial state:")
        
        for init_name, init_state in initial_states.items():
            final_state = chosen_gate_matrix @ init_state
            # This printout fulfills the "output each number in the final equation" requirement.
            print(f"{chosen_gate_name} {init_name}  ->  {state_to_string(final_state)}")
            
        print("\nNone of the resulting states are |i> or |-i>, so no deaths occur.")
    else:
        print("No safe gate found among the options.")

solve_trolley_problem()
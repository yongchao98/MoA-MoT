import numpy as np

def are_states_equal(state1, state2, tolerance=1e-9):
    """Checks if two state vectors are physically equivalent (proportional)."""
    # The absolute value of the inner product of two normalized vectors must be 1.
    norm1 = np.linalg.norm(state1)
    norm2 = np.linalg.norm(state2)
    if norm1 < tolerance or norm2 < tolerance:
        return abs(norm1 - norm2) < tolerance

    s1_normalized = state1 / norm1
    s2_normalized = state2 / norm2
    
    inner_product_mag = abs(np.vdot(s1_normalized, s2_normalized))
    return np.isclose(inner_product_mag, 1.0, atol=tolerance)

def run_analysis():
    """
    Analyzes the quantum trolley problem to find the safe gate operation.
    """
    # Define basis states using the computational basis |0> and |1>
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)
    
    # The six standard basis states
    s_minus = (s0 - s1) / np.sqrt(2)
    s_i = (s0 + 1j * s1) / np.sqrt(2)
    s_minus_i = (s0 - 1j * s1) / np.sqrt(2)

    # The problem specifies the initial state is one of these five
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": s_minus,
        "|i>": s_i,
        "|-i>": s_minus_i,
    }

    # The two fatal outcomes (left and right tracks)
    death_states = {"|i>": s_i, "|-i>": s_minus_i}

    # A selection of well-defined, single-qubit gates from the choice list
    gates = {
        "B. T†": np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]]),
        "C. S†": np.array([[1, 0], [0, -1j]]),
        "H. Z": np.array([[1, 0], [0, -1]]),
        "L. X": np.array([[0, 1], [1, 0]]),
        "M. Y": np.array([[0, -1j], [1j, 0]]),
        "T. H": 1/np.sqrt(2) * np.array([[1, 1], [1, -1]]),
        "U. I": np.identity(2, dtype=complex),
        "W. S": np.array([[1, 0], [0, 1j]]),
        "Y. T": np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]]),
    }

    print("--- Analysis of Quantum Gate Operations ---")
    solutions = []
    for name, gate_matrix in gates.items():
        is_safe = True
        failure_reason = ""
        for init_name, init_state in initial_states.items():
            final_state = gate_matrix @ init_state
            for death_name, death_state in death_states.items():
                if are_states_equal(final_state, death_state):
                    is_safe = False
                    failure_reason = f"Gate {name} fails because it maps initial state {init_name} to fatal state {death_name}."
                    break
            if not is_safe:
                break
        
        if is_safe:
            print(f"SAFE: {name}")
            solutions.append(name.split('.')[0])
        else:
            print(f"FAIL: {failure_reason}")

    print("\n--- Equation for a Failing Gate (C. S†) ---")
    # S† |−⟩ = S† * (1/√2 * (|0⟩ - |1⟩))
    #         = 1/√2 * (S†|0⟩ - S†|1⟩)
    #         = 1/√2 * (|0⟩ - (-i)|1⟩)
    #         = 1/√2 * (|0⟩ + i|1⟩)
    #         = |i⟩
    print("Equation: S† |−⟩ = |i⟩")
    print("In this equation, the numbers are:")
    # Coefficients: 1/sqrt(2), 1, 1, 1/sqrt(2), 1, 1, 1/sqrt(2), -i, 1, 1/sqrt(2), 1, i, 1
    # For simplicity, printing the non-trivial components and values.
    # The coefficients of the final state are 1/sqrt(2) and i/sqrt(2).
    print("Initial state |−⟩ is (1/√2)|0⟩ + (-1/√2)|1⟩")
    print("S† acts on |1⟩, changing its phase by -i")
    print("The final state is (1/√2)|0⟩ + (i/√2)|1⟩, which is the definition of |i⟩.")
    print("\nAs the simulation and the equation show, applying S† when the lever is in the |−⟩ state leads to the |i⟩ state (right track), causing deaths.")
    print(f"\nThe only operations that are safe for all possible initial states are: {', '.join(solutions)}")


run_analysis()
<<<Y>>>
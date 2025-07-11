import numpy as np

def ket_to_string(ket, precision=3):
    """Converts a ket vector to a readable string like 'a|0> + b|1>'."""
    ket = np.asarray(ket).flatten()
    val0 = ket[0]
    val1 = ket[1]
    
    # Handle the first term
    if np.isclose(val0, 0):
        term0_str = ""
    else:
        term0_str = f"({val0:.{precision}f})|0⟩"

    # Handle the second term
    if np.isclose(val1, 0):
        term1_str = ""
    else:
        if val1.imag < 0:
            op = "-"
            val1_abs = -val1
        else:
            op = "+"
            val1_abs = val1

        # Don't print the plus sign if the first term is zero
        if term0_str == "":
             term1_str = f"({val1:.{precision}f})|1⟩"
        else:
             term1_str = f" {op} ({val1_abs:.{precision}f})|1⟩"

    if not term0_str and not term1_str:
        return "0"
        
    return term0_str + term1_str

def is_proportional(v1, v2, tol=1e-9):
    """Checks if vector v1 is proportional to vector v2."""
    v1 = np.asarray(v1, dtype=complex).flatten()
    v2 = np.asarray(v2, dtype=complex).flatten()
    # Ensure they are normalized
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    # Check if the absolute value of their inner product is close to 1
    dot_product = np.vdot(v1, v2)
    return np.isclose(np.abs(dot_product), 1.0, atol=tol)

def solve_trolley_problem():
    """
    Solves the quantum trolley problem by finding a safe gate operation.
    """
    # --- 1. Define the states ---
    sqrt2 = np.sqrt(2)
    
    # Initial states (ket vectors)
    ket_0 = np.array([1, 0], dtype=complex)
    ket_1 = np.array([0, 1], dtype=complex)
    ket_minus = np.array([1/sqrt2, -1/sqrt2], dtype=complex)
    ket_i = np.array([1/sqrt2, 1j/sqrt2], dtype=complex)
    ket_minus_i = np.array([1/sqrt2, -1j/sqrt2], dtype=complex)

    initial_states = {
        "|0⟩": ket_0,
        "|1⟩": ket_1,
        "|-⟩": ket_minus,
        "|i⟩": ket_i,
        "|-i⟩": ket_minus_i,
    }
    
    # "Death" states
    death_states = [ket_i, ket_minus_i]

    # --- 2. Define and select the gate ---
    # The T gate matrix
    # T represents a rotation around the Z-axis by pi/4
    t_gate = np.array([
        [1, 0],
        [0, np.exp(1j * np.pi / 4)]
    ], dtype=complex)
    
    chosen_gate = t_gate
    chosen_gate_name = "T"
    
    print(f"Chosen Action: Apply the '{chosen_gate_name}' gate.")
    print("Verifying this action for all possible initial states...")
    print("-" * 50)
    
    is_safe_gate = True
    
    # --- 3. Apply the gate and check outcomes ---
    for name, initial_ket in initial_states.items():
        # Apply the gate: final_ket = U * initial_ket
        final_ket = chosen_gate @ initial_ket
        
        # Check if the final state is a death state
        is_death = False
        for death_state in death_states:
            if is_proportional(final_ket, death_state):
                is_death = True
                break
        
        print(f"Initial State: {name}")
        print(f"Result of T{name}: {ket_to_string(final_ket)}")
        
        if is_death:
            print("Outcome: UNSAFE. This state leads to a death state.")
            is_safe_gate = False
        else:
            print("Outcome: SAFE. This state does not lead to a death state.")
        print()

    print("-" * 50)
    
    if is_safe_gate:
        print("Conclusion: The 'T' gate is a safe operation that avoids all human deaths.")
    else:
        print("Conclusion: The 'T' gate is NOT a safe operation.")


solve_trolley_problem()
<<<Y>>>
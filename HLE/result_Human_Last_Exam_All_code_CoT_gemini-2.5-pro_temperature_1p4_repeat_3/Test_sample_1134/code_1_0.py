import numpy as np

def is_proportional(v1, v2, tolerance=1e-9):
    """Checks if two state vectors are proportional (represent the same state)."""
    # Normalize both vectors
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 < tolerance or norm2 < tolerance:
        return False
    v1_norm = v1 / norm1
    v2_norm = v2 / norm2
    # Check if the absolute value of their inner product is close to 1
    return abs(np.abs(np.vdot(v1_norm, v2_norm)) - 1.0) < tolerance

def state_to_string(v, tolerance=1e-9):
    """Converts a state vector to its string representation if it's a known basis state."""
    v_map = {
        "|0>": np.array([1, 0]),
        "|1>": np.array([0, 1]),
        "|+>": (1/np.sqrt(2)) * np.array([1, 1]),
        "|->": (1/np.sqrt(2)) * np.array([1, -1]),
        "|i>": (1/np.sqrt(2)) * np.array([1, 1j]),
        "|-i>": (1/np.sqrt(2)) * np.array([1, -1j])
    }
    for name, basis_v in v_map.items():
        if is_proportional(v, basis_v, tolerance):
            return name
    return f"[{v[0]:.2f}, {v[1]:.2f}]"

def main():
    # Define basis states
    states = {
        "|0>": np.array([1, 0], dtype=complex),
        "|1>": np.array([0, 1], dtype=complex),
        "|->": (1/np.sqrt(2)) * np.array([1, -1], dtype=complex),
        "|i>": (1/np.sqrt(2)) * np.array([1, 1j], dtype=complex),
        "|-i>": (1/np.sqrt(2)) * np.array([1, -1j], dtype=complex),
    }

    # States that must be avoided
    death_states = {
        "|i>": states["|i>"],
        "|-i>": states["|-i>"],
    }

    # Define quantum gates from the list
    gates = {
        "H": (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
        "S": np.array([[1, 0], [0, 1j]], dtype=complex),
        "S†": np.array([[1, 0], [0, -1j]], dtype=complex),
        "T": np.array([[1, 0], [0, np.exp(1j*np.pi/4)]], dtype=complex),
        "T†": np.array([[1, 0], [0, np.exp(-1j*np.pi/4)]], dtype=complex),
        "I": np.identity(2, dtype=complex),
    }

    safe_gates = []

    print("--- Analyzing Quantum Gate Safety ---\n")

    for gate_name, gate_matrix in gates.items():
        is_gate_safe = True
        print(f"Testing Gate: {gate_name}")
        print("-" * 25)
        for state_name, state_vector in states.items():
            # Apply the gate
            result_vector = gate_matrix @ state_vector

            # Check if the result is a death state
            is_fatal = False
            for death_name, death_vector in death_states.items():
                if is_proportional(result_vector, death_vector):
                    is_fatal = True
                    break
            
            # Print the equation and the result
            final_state_str = state_to_string(result_vector)
            print(f"Equation: {gate_name} {state_name} = {final_state_str}")
            if final_state_str.startswith("["): # If it's not a named state, show vector
                 print(f"Result Vector: [{result_vector[0]:.3f} {result_vector[1]:.3f}]")

            if is_fatal:
                print(f"Result: DANGER! This leads to the {final_state_str} state.\n")
                is_gate_safe = False
                break # No need to check other initial states for this gate
            else:
                print("Result: Safe.\n")

        if is_gate_safe:
            safe_gates.append(gate_name)
        print("="*40 + "\n")

    print("\n--- Conclusion ---")
    if safe_gates:
        print(f"The following gates are safe for all possible initial states: {', '.join(safe_gates)}")
        print("Applying any of these gates will guarantee no one is harmed.")
    else:
        print("No safe gate was found in the analysis.")

if __name__ == "__main__":
    main()
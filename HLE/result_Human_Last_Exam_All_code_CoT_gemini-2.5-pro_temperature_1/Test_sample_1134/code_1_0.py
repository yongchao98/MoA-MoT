import numpy as np

def normalize(v):
    """Normalizes a vector."""
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def are_equivalent(v1, v2, tolerance=1e-9):
    """Checks if two state vectors are physically equivalent (equal up to a global phase)."""
    v1 = normalize(np.array(v1, dtype=complex))
    v2 = normalize(np.array(v2, dtype=complex))
    
    # Find the first non-zero element to determine the phase difference
    phase_v1 = 0
    phase_v2 = 0
    for i in range(len(v1)):
        if abs(v1[i]) > tolerance:
            phase_v1 = v1[i]
            phase_v2 = v2[i]
            break
            
    if abs(phase_v1) < tolerance: # v1 is a zero vector
        return np.allclose(v1, v2)

    # Calculate phase factor c such that v2 = c * v1
    # We expect phase_v2 = c * phase_v1 => c = phase_v2 / phase_v1
    c = phase_v2 / phase_v1
    
    # Check if magnitude of c is 1 and if v2 is indeed c * v1
    return abs(abs(c) - 1) < tolerance and np.allclose(v2, c * v1)

def state_to_string(v):
    """Converts a state vector to its string representation."""
    v = normalize(v)
    # Check against basis states
    if are_equivalent(v, S["|0>"]): return "|0>"
    if are_equivalent(v, S["|1>"]): return "|1>"
    if are_equivalent(v, S["|+>"]): return "|+>"
    if are_equivalent(v, S["|->"]): return "|->"
    if are_equivalent(v, S["|i>"]): return "|i>"
    if are_equivalent(v, S["|-i>"]): return "|-i>"
    return f"[{v[0]:.2f}, {v[1]:.2f}]"

def get_phase_string(v_final, v_original):
    """Gets the string for the global phase c in v_final = c * v_original."""
    v_final = normalize(np.array(v_final, dtype=complex))
    v_original = normalize(np.array(v_original, dtype=complex))
    # Find first non-zero element to calculate phase
    phase_final = 0
    phase_original = 0
    for i in range(len(v_final)):
        if abs(v_final[i]) > 1e-9:
            phase_final = v_final[i]
            phase_original = v_original[i]
            break
    c = phase_final / phase_original
    if np.allclose(c, 1):
        return ""
    return f"({c:.2f}) "


# --- 1. Define States ---
S = {
    "|0>": np.array([1, 0], dtype=complex),
    "|1>": np.array([0, 1], dtype=complex),
    "|+>": 1/np.sqrt(2) * np.array([1, 1], dtype=complex),
    "|->": 1/np.sqrt(2) * np.array([1, -1], dtype=complex),
    "|i>": 1/np.sqrt(2) * np.array([1, 1j], dtype=complex),
    "|-i>": 1/np.sqrt(2) * np.array([1, -1j], dtype=complex),
}

initial_states = {
    "|0>": S["|0>"],
    "|1>": S["|1>"],
    "|->": S["|->"],
    "|i>": S["|i>"],
    "|-i>": S["|-i>"],
}

death_states = {
    "|i>": S["|i>"],
    "|-i>": S["|-i>"],
}

# --- 2. Define Gates ---
GATES = {
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "H": 1/np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=complex),
    "S": np.array([[1, 0], [0, 1j]], dtype=complex),
    "S†": np.array([[1, 0], [0, -1j]], dtype=complex),
    "T": np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex),
    "T†": np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=complex),
}

# --- 3. Test Gates ---
solution = None
print("Analyzing potential operations...\n")

for gate_name, gate_matrix in GATES.items():
    is_safe = True
    print(f"--- Testing Gate: {gate_name} ---")
    for init_name, init_vec in initial_states.items():
        final_vec = gate_matrix @ init_vec
        
        crashed = False
        for death_name, death_vec in death_states.items():
            if are_equivalent(final_vec, death_vec):
                phase_str = get_phase_string(final_vec, death_vec)
                print(f"FAILED: Applying {gate_name} to {init_name} results in death.")
                print(f"   Equation: {gate_name}{init_name} = {phase_str}{death_name}")
                is_safe = False
                crashed = True
                break # No need to check other death state
        if crashed:
            break # No need to check other initial states
    
    if is_safe:
        solution = gate_name
        print(f"SUCCESS: Gate {gate_name} is safe for all initial states.")

    print("-" * (20 + len(gate_name)))
    print()

if solution:
    print(f"\nConclusion: The only safe operation is {solution}.")
else:
    print("\nConclusion: No single safe operation found among the tested gates.")

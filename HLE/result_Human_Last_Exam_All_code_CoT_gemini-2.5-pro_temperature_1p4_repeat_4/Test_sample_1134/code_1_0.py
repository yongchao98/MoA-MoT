import numpy as np

def format_state_vector(v):
    """Formats a 2x1 complex vector for printing."""
    # Normalize for consistent representation, although T-dagger is unitary
    norm = np.linalg.norm(v)
    if norm > 1e-9:
        v = v / norm
        
    # Find a global phase to make the first non-zero element real and positive
    phase = 1.0
    for x in v:
        if not np.isclose(x, 0):
            phase = np.exp(-1j * np.angle(x))
            if np.isclose(x * phase, -1.0):
                 phase = -phase
            break
    v_phased = v * phase
    
    # Clean up small numbers to zero
    v_phased[np.isclose(v_phased, 0)] = 0
    
    # Format components
    c0_str = f"{v_phased[0]:.3f}".replace('+0.000j', '').replace('-0.000j', '').replace('j', 'i')
    c1_str = f"{v_phased[1]:.3f}".replace('+0.000j', '').replace('-0.000j', '').replace('j', 'i')
    
    return f"[{c0_str}, {c1_str}]ᵀ"


def is_death_state(state_vector):
    """Checks if a state vector is |i> or |-i> up to a global phase."""
    v = state_vector / np.linalg.norm(state_vector)
    i_state = (1/np.sqrt(2)) * np.array([1, 1j])
    neg_i_state = (1/np.sqrt(2)) * np.array([1, -1j])
    
    is_i = np.allclose(np.abs(np.dot(i_state.conj(), v)), 1.0)
    is_neg_i = np.allclose(np.abs(np.dot(neg_i_state.conj(), v)), 1.0)
    
    return is_i or is_neg_i

# Define the quantum states and gate T_dagger (T†)
states = {
    "|0>": np.array([1, 0], dtype=np.complex128),
    "|1>": np.array([0, 1], dtype=np.complex128),
    "|->": (1/np.sqrt(2)) * np.array([1, -1], dtype=np.complex128),
    "|i>": (1/np.sqrt(2)) * np.array([1, 1j], dtype=np.complex128),
    "|-i>": (1/np.sqrt(2)) * np.array([1, -1j], dtype=np.complex128)
}

# Gate T-dagger (B)
T_dagger = np.array([[1, 0], [0, np.exp(-1j * np.pi/4)]], dtype=np.complex128)
gate_name = "T†"

# Main logic to apply the gate and print equations
print(f"Applying the {gate_name} gate to all possible initial states:\n")

all_safe = True
for name, vector in states.items():
    result_vector = T_dagger @ vector
    
    # Format initial state vector for printing
    initial_vec_str = format_state_vector(vector)
    # Format final state vector for printing
    final_vec_str = format_state_vector(result_vector)

    print(f"{gate_name} {name} = {gate_name} * {initial_vec_str} = {final_vec_str}")
    
    if is_death_state(result_vector):
        print("--> This is a DEATH state. The action is not safe.\n")
        all_safe = False
    else:
        print("--> This is a SAFE state.\n")

if all_safe:
    print("Conclusion: Applying the T† gate successfully avoids all human deaths.")

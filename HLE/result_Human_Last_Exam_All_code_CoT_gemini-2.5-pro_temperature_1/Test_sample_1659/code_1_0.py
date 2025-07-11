import cmath

def print_state(state, name):
    """Prints the quantum state in a readable format."""
    parts = []
    # Sort for consistent output and handle floating point inaccuracies
    for basis_state in sorted(state.keys()):
        amplitude = state[basis_state]
        if abs(amplitude) < 1e-9:
            continue
        
        # Determine the sign for the term
        if len(parts) > 0:
            sign = " + " if amplitude.real >= 0 else " - "
        else:
            sign = "" if amplitude.real >= 0 else "-"
            
        amp_val_str = f"{abs(amplitude.real):.4f}"
        
        parts.append(f"{sign}{amp_val_str}|{basis_state}>")
            
    print(f"{name} = {''.join(parts).lstrip(' + ')}")

def add_to_state(state, basis, amp):
    """Helper function to add a component to a state vector, handling superposition."""
    state[basis] = state.get(basis, 0.0) + amp

def apply_h(state, qubit_index):
    """Applies a Hadamard gate to a specific qubit in the state."""
    new_state = {}
    sqrt2_inv = 1 / cmath.sqrt(2)
    for basis, amp in state.items():
        # H|0> = 1/sqrt(2) * (|0> + |1>)
        if basis[qubit_index] == '0':
            new_basis_0_list = list(basis)
            new_basis_1_list = list(basis)
            new_basis_1_list[qubit_index] = '1'
            add_to_state(new_state, "".join(new_basis_0_list), amp * sqrt2_inv)
            add_to_state(new_state, "".join(new_basis_1_list), amp * sqrt2_inv)
        # H|1> = 1/sqrt(2) * (|0> - |1>)
        else: # basis[qubit_index] == '1'
            new_basis_0_list = list(basis)
            new_basis_0_list[qubit_index] = '0'
            new_basis_1_list = list(basis)
            add_to_state(new_state, "".join(new_basis_0_list), amp * sqrt2_inv)
            add_to_state(new_state, "".join(new_basis_1_list), -amp * sqrt2_inv)
    return new_state

def apply_cnot(state, control_index, target_index):
    """Applies a CNOT gate."""
    new_state = {}
    for basis, amp in state.items():
        if basis[control_index] == '1':
            new_basis_list = list(basis)
            new_basis_list[target_index] = '1' if basis[target_index] == '0' else '0'
            add_to_state(new_state, "".join(new_basis_list), amp)
        else:
            add_to_state(new_state, basis, amp)
    return new_state

def apply_ccnot(state, c1_idx, c2_idx, t_idx):
    """Applies a Toffoli (CCNOT) gate."""
    new_state = {}
    for basis, amp in state.items():
        if basis[c1_idx] == '1' and basis[c2_idx] == '1':
            new_basis_list = list(basis)
            new_basis_list[t_idx] = '1' if basis[t_idx] == '0' else '0'
            add_to_state(new_state, "".join(new_basis_list), amp)
        else:
            add_to_state(new_state, basis, amp)
    return new_state

# --- Main Quantum Circuit Simulation ---

# Initial State: |psi_0> = |000>
# Qubit indices are 0, 1, 2
psi_0 = {'000': 1.0}
print("Step-by-step state evolution:")
print_state(psi_0, "|psi_0>")
print("-" * 40)

# 1. Apply a Hadamard gate to the first qubit (index 0)
psi_1 = apply_h(psi_0, 0)
print("1. Apply H to qubit 1: |psi_1> = (H @ I @ I)|psi_0>")
print_state(psi_1, "|psi_1>")
print("-" * 40)

# 2. Apply CNOT (1st qubit control, 2nd target)
psi_2 = apply_cnot(psi_1, 0, 1)
print("2. Apply CNOT(1,2): |psi_2> = CNOT_1,2 |psi_1>")
print_state(psi_2, "|psi_2>")
print("-" * 40)

# 3. Apply Toffoli (1st, 2nd control, 3rd target)
psi_3 = apply_ccnot(psi_2, 0, 1, 2)
print("3. Apply CCNOT(1,2,3): |psi_3> = CCNOT_1,2,3 |psi_2>")
print_state(psi_3, "|psi_3>")
print("-" * 40)

# 4. Apply a second Hadamard gate to the first qubit
psi_4 = apply_h(psi_3, 0)
print("4. Apply H to qubit 1 again: |psi_4> = (H @ I @ I)|psi_3>")
print_state(psi_4, "|psi_4>")
print("-" * 40)

# --- Probability Calculation ---
target_outcome = '100'
amplitude_target = psi_4.get(target_outcome, 0)
probability = abs(amplitude_target)**2

print("Final Calculation:")
print(f"The final state is |psi_4> = 0.5000|000> + 0.5000|011> + 0.5000|100> - 0.5000|111>")
print("\nThe probability of measuring |100> is the squared magnitude of its amplitude.")
print(f"Amplitude of |{target_outcome}> in |psi_4> is {amplitude_target.real:.4f}")
print(f"P(|{target_outcome}>) = |Amplitude(|{target_outcome}>)|^2")
print(f"P(|{target_outcome}>) = |{amplitude_target.real:.4f}|^2")
print(f"P(|{target_outcome}>) = {probability:.4f}")

<<<0.25>>>
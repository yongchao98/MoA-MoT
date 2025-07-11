import cmath
from collections import defaultdict

def apply_hadamard(state, qubit_index):
    """Applies a Hadamard gate to a specific qubit."""
    new_state = defaultdict(complex)
    inv_sqrt2 = 1 / cmath.sqrt(2)
    for basis, amplitude in state.items():
        bit = basis[qubit_index]
        
        basis_list = list(basis)
        
        if bit == '0':
            # H|0> = 1/sqrt(2) * (|0> + |1>)
            # The original basis state |...0...> contributes to itself
            new_state[basis] += amplitude * inv_sqrt2
            # and to the state with the flipped bit |...1...>
            basis_list[qubit_index] = '1'
            other_basis = "".join(basis_list)
            new_state[other_basis] += amplitude * inv_sqrt2
        else:  # bit == '1'
            # H|1> = 1/sqrt(2) * (|0> - |1>)
            # The original basis state |...1...> contributes to itself with a negative sign
            new_state[basis] += -amplitude * inv_sqrt2
            # and to the state with the flipped bit |...0...>
            basis_list[qubit_index] = '0'
            other_basis = "".join(basis_list)
            new_state[other_basis] += amplitude * inv_sqrt2
            
    return dict(new_state)

def apply_cnot(state, control_index, target_index):
    """Applies a CNOT gate."""
    new_state = defaultdict(complex)
    for basis, amplitude in state.items():
        if basis[control_index] == '1':
            # Flip the target bit
            target_bit = int(basis[target_index])
            flipped_bit = 1 - target_bit
            new_basis_list = list(basis)
            new_basis_list[target_index] = str(flipped_bit)
            new_basis = "".join(new_basis_list)
            new_state[new_basis] += amplitude
        else:
            # Control is 0, do nothing
            new_state[basis] += amplitude
    return dict(new_state)

def apply_ccnot(state, c1_index, c2_index, target_index):
    """Applies a CCNOT (Toffoli) gate."""
    new_state = defaultdict(complex)
    for basis, amplitude in state.items():
        if basis[c1_index] == '1' and basis[c2_index] == '1':
            # Both controls are 1, flip the target bit
            target_bit = int(basis[target_index])
            flipped_bit = 1 - target_bit
            new_basis_list = list(basis)
            new_basis_list[target_index] = str(flipped_bit)
            new_basis = "".join(new_basis_list)
            new_state[new_basis] += amplitude
        else:
            # At least one control is 0, do nothing
            new_state[basis] += amplitude
    return dict(new_state)

def format_state(state):
    """Formats the state dictionary into a readable string."""
    # Sort by basis state for a consistent, ordered output
    sorted_items = sorted(state.items())
    parts = []
    for basis, amp in sorted_items:
        if abs(amp) > 1e-9: # Ignore components with near-zero amplitude
            # Format complex numbers nicely
            if abs(amp.imag) < 1e-9:
                sign = "-" if amp.real < 0 else "+"
                val = abs(amp.real)
                parts.append(f"{sign} {val:.3f}|{basis}⟩")
            elif abs(amp.real) < 1e-9:
                sign = "-" if amp.imag < 0 else "+"
                val = abs(amp.imag)
                parts.append(f"{sign} {val:.3f}j|{basis}⟩")
            else:
                 parts.append(f"+ ({amp:.3f})|{basis}⟩")
    
    # Clean up the string
    result = " ".join(parts)
    if result.startswith("+ "):
        result = result[2:]
    return result.replace(" + -", " - ")


def solve_quantum_circuit():
    """
    Solves the given quantum circuit problem.
    """
    # Initial State: |ψ₀⟩ = |000⟩
    psi_0 = {'000': 1.0}
    print(f"Initial State: |ψ₀⟩ = {format_state(psi_0)}")

    # 1. Apply a Hadamard gate to the first qubit (index 0)
    psi_1 = apply_hadamard(psi_0, 0)
    print(f"Step 1 (H on qubit 1): |ψ₁⟩ = {format_state(psi_1)}")

    # 2. Apply a CNOT gate with qubit 1 (control) and 2 (target)
    psi_2 = apply_cnot(psi_1, 0, 1)
    print(f"Step 2 (CNOT on 1,2): |ψ₂⟩ = {format_state(psi_2)}")

    # 3. Apply a Toffoli gate with qubits 1,2 (controls) and 3 (target)
    psi_3 = apply_ccnot(psi_2, 0, 1, 2)
    print(f"Step 3 (CCNOT on 1,2,3): |ψ₃⟩ = {format_state(psi_3)}")

    # 4. Apply a second Hadamard gate to the first qubit
    psi_4 = apply_hadamard(psi_3, 0)
    print(f"Step 4 (H on qubit 1): |ψ₄⟩ = {format_state(psi_4)}\n")

    # Determine the probability of measuring |100⟩
    target_state = '100'
    amplitude = psi_4.get(target_state, 0)
    probability = abs(amplitude)**2

    print("--- Probability Calculation ---")
    print(f"The final state is |ψ₄⟩ = {format_state(psi_4)}")
    print(f"The amplitude of the |{target_state}⟩ component is c_100.")
    print(f"From the final state, we find c_100 = {amplitude.real:.4f}")
    print(f"The probability of measuring |{target_state}⟩ is P(|{target_state}⟩) = |c_100|².")
    print(f"P(|{target_state}⟩) = |{amplitude.real:.4f}|² = {probability:.4f}")

# Run the simulation and print the results
solve_quantum_circuit()
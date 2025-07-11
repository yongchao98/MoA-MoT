import math
from collections import defaultdict

def format_coeff(coeff, use_unicode=True):
    """Formats a single coefficient into a nice string (e.g., 1/2, 1/√2)."""
    val = abs(coeff)
    sign = "-" if coeff < 0 else ""
    if abs(val - 1.0) < 1e-9: return f"{sign}1"
    if abs(val - 0.5) < 1e-9: return f"{sign}1/2"
    if abs(val - 0.25) < 1e-9: return f"{sign}1/4"
    # Python's math.sqrt(2) might not be exactly what the unicode char represents
    # but it's close enough for display.
    if abs(val - 1/math.sqrt(2)) < 1e-9: return f"{sign}1/{'√2' if use_unicode else 'sqrt(2)'}"
    return f"{coeff:.4f}"

def format_state(state, use_unicode=True):
    """Formats a state dictionary into a human-readable ket notation string."""
    if not state:
        return "0"
    parts = []
    # Sort by ket for consistent output (e.g., '000', '001', ...)
    for ket in sorted(state.keys()):
        coeff = state[ket]
        if abs(coeff) < 1e-9:
            continue
        
        sign = " - " if coeff < 0 else " + "
        val = abs(coeff)
        
        # Format the coefficient's absolute value
        if abs(val - 1.0) < 1e-9:
            coeff_str = ""
        else:
            # Add parenthesis around the coefficient
            coeff_str = f"({format_coeff(val, use_unicode)})"
        
        term = f"{coeff_str}|{ket}⟩"
        
        # For the very first term, handle the sign separately
        if not parts:
            if coeff < 0:
                parts.append(f"-{term}")
            else:
                parts.append(term)
        else:
            parts.append(sign)
            parts.append(term)
    
    return "".join(parts)

def apply_H(state, qubit_index):
    """Applies a Hadamard gate to a specific qubit in the state."""
    new_state = defaultdict(float)
    sqrt2_inv = 1 / math.sqrt(2)
    for ket, coeff in state.items():
        bit = ket[qubit_index]
        ket_list = list(ket)
        if bit == '0':
            # |0> -> 1/√2 * (|0> + |1>)
            new_state[ket] += coeff * sqrt2_inv
            ket_list[qubit_index] = '1'
            new_ket = "".join(ket_list)
            new_state[new_ket] += coeff * sqrt2_inv
        else: # bit == '1'
            # |1> -> 1/√2 * (|0> - |1>)
            ket_list[qubit_index] = '0'
            new_ket = "".join(ket_list)
            new_state[new_ket] += coeff * sqrt2_inv
            new_state[ket] += -coeff * sqrt2_inv
    return dict(new_state)

def apply_CNOT(state, control_index, target_index):
    """Applies a CNOT gate to the state."""
    new_state = defaultdict(float)
    for ket, coeff in state.items():
        if ket[control_index] == '1':
            ket_list = list(ket)
            target_bit = ket_list[target_index]
            ket_list[target_index] = '1' if target_bit == '0' else '0'
            new_ket = "".join(ket_list)
            new_state[new_ket] += coeff
        else:
            new_state[ket] += coeff
    return dict(new_state)

def apply_CCNOT(state, control1_index, control2_index, target_index):
    """Applies a CCNOT (Toffoli) gate to the state."""
    new_state = defaultdict(float)
    for ket, coeff in state.items():
        if ket[control1_index] == '1' and ket[control2_index] == '1':
            ket_list = list(ket)
            target_bit = ket_list[target_index]
            ket_list[target_index] = '1' if target_bit == '0' else '0'
            new_ket = "".join(ket_list)
            new_state[new_ket] += coeff
        else:
            new_state[ket] += coeff
    return dict(new_state)

def solve_quantum_circuit():
    """
    Solves the given quantum circuit problem, printing each step and the final probability.
    Note: Qubit numbering in print statements (1,2,3) matches the problem description,
          while the underlying code uses 0-based indexing (0,1,2).
    """
    print("Step-by-step calculation of the quantum state:")
    print("-" * 60)

    # Step 0: Initial State
    psi_0 = {'000': 1.0}
    print(f"Initial State: |ψ₀⟩ = {format_state(psi_0)}\n")

    # Step 1: Apply Hadamard to the first qubit (index 0)
    psi_1 = apply_H(psi_0, 0)
    print(f"Step 1: Apply H to qubit 1 -> |ψ₁⟩ = (H ⊗ I ⊗ I)|ψ₀⟩")
    print(f"|ψ₁⟩ = {format_state(psi_1)}\n")

    # Step 2: Apply CNOT with qubit 1 (index 0) as control and qubit 2 (index 1) as target
    psi_2 = apply_CNOT(psi_1, 0, 1)
    print(f"Step 2: Apply CNOT(1,2) -> |ψ₂⟩ = CNOT₁,₂|ψ₁⟩")
    print(f"|ψ₂⟩ = {format_state(psi_2)}\n")

    # Step 3: Apply Toffoli gate with qubits 1,2 as control and qubit 3 as target
    psi_3 = apply_CCNOT(psi_2, 0, 1, 2)
    print(f"Step 3: Apply CCNOT(1,2,3) -> |ψ₃⟩ = CCNOT₁,₂,₃|ψ₂⟩")
    print(f"|ψ₃⟩ = {format_state(psi_3)}\n")

    # Step 4: Apply Hadamard to the first qubit (index 0)
    psi_4 = apply_H(psi_3, 0)
    print(f"Step 4: Apply H to qubit 1 -> |ψ₄⟩ = (H ⊗ I ⊗ I)|ψ₃⟩")
    print(f"|ψ₄⟩ = {format_state(psi_4)}\n")

    print("-" * 60)
    print("Final Probability Calculation:")

    target_ket = '100'
    amplitude = psi_4.get(target_ket, 0)
    probability = abs(amplitude)**2

    print(f"The final state is |ψ₄⟩ = {format_state(psi_4)}.")
    print(f"The amplitude (coefficient) of the target state |{target_ket}⟩ is {format_coeff(amplitude)}.")
    print(f"The probability of measuring |{target_ket}⟩ is the square of the magnitude of its amplitude.")
    print(f"P(|{target_ket}⟩) = |Amplitude(|{target_ket}⟩)|² = |{format_coeff(amplitude)}|² = ({amplitude:.2f}) * ({amplitude:.2f}) = {probability}")

# Execute the simulation
solve_quantum_circuit()
<<<0.25>>>
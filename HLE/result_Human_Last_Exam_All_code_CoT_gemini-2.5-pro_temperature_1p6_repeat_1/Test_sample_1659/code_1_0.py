import math
from collections import defaultdict

def apply_gate(state, gate_function, *args):
    """
    Applies a quantum gate to a state by iterating over its basis states.

    Args:
        state (dict): The quantum state as {basis_state_str: amplitude}.
        gate_function (function): The function implementing the gate's logic.
        *args: Qubit indices to be passed to the gate_function.

    Returns:
        dict: The new quantum state after applying the gate.
    """
    new_state = defaultdict(complex)
    for basis_state, amplitude in state.items():
        if abs(amplitude) < 1e-9:  # Ignore zero-amplitude states
            continue
        # A gate maps one basis state to a superposition of basis states.
        # The gate function returns a list of: (new_basis_state, amplitude_multiplier)
        outcomes = gate_function(basis_state, *args)
        for new_basis_state, amp_multiplier in outcomes:
            new_state[new_basis_state] += amplitude * amp_multiplier
    # Filter out any resulting zero-amplitude states
    return {k: v for k, v in new_state.items() if abs(v) > 1e-9}

def hadamard_on_basis(basis_state, qubit_idx):
    """
    Simulates applying a Hadamard gate to a qubit in a basis state.
    Returns a list of [(new_state, amplitude_multiplier), ...].
    """
    qubit_val = int(basis_state[qubit_idx])
    prefix = basis_state[:qubit_idx]
    suffix = basis_state[qubit_idx + 1:]
    
    amp_multiplier = 1 / math.sqrt(2)
    
    state0_str = prefix + '0' + suffix
    state1_str = prefix + '1' + suffix
    
    if qubit_val == 0: # |0> -> 1/sqrt(2) * (|0> + |1>)
        return [(state0_str, amp_multiplier), (state1_str, amp_multiplier)]
    else: # |1> -> 1/sqrt(2) * (|0> - |1>)
        return [(state0_str, amp_multiplier), (state1_str, -amp_multiplier)]

def cnot_on_basis(basis_state, control_idx, target_idx):
    """
    Simulates applying a CNOT gate to two qubits in a basis state.
    """
    if basis_state[control_idx] == '1':
        # Flip the target qubit
        state_list = list(basis_state)
        target_val = int(state_list[target_idx])
        state_list[target_idx] = str(1 - target_val)
        new_basis_state = "".join(state_list)
        return [(new_basis_state, 1.0)]
    else:
        # Control is 0, state is unchanged
        return [(basis_state, 1.0)]

def toffoli_on_basis(basis_state, c1_idx, c2_idx, target_idx):
    """
    Simulates applying a Toffoli (CCNOT) gate to three qubits in a basis state.
    """
    if basis_state[c1_idx] == '1' and basis_state[c2_idx] == '1':
        # Flip the target qubit
        state_list = list(basis_state)
        target_val = int(state_list[target_idx])
        state_list[target_idx] = str(1 - target_val)
        new_basis_state = "".join(state_list)
        return [(new_basis_state, 1.0)]
    else:
        # Controls are not both 1, state is unchanged
        return [(basis_state, 1.0)]

def format_state(state):
    """Formats a quantum state dictionary into a readable string."""
    parts = []
    # Sort by basis state string for consistent output
    sorted_items = sorted(state.items())
    
    for i, (basis, amp) in enumerate(sorted_items):
        if abs(amp) < 1e-9: continue
        
        sign = " - " if amp.real < 0 else " + "
        if i == 0 and amp.real > 0: sign = ""
        elif i == 0 and amp.real < 0: sign = "-"

        amp_val = abs(amp.real)
        amp_str = f"{amp_val:.4f}"
        
        parts.append(f"{sign}{amp_str}|{basis}>")
        
    return "".join(parts)


if __name__ == '__main__':
    # Initial State: |psi_0> = |000>
    psi_0 = {'000': 1.0}
    print(f"Step 0: Initial state |psi_0> = {format_state(psi_0)}")

    # Step 1: Apply H to the first qubit (index 0)
    # |psi_1> = (H @ I @ I) |psi_0>
    psi_1 = apply_gate(psi_0, hadamard_on_basis, 0)
    print(f"Step 1: After H on qubit 0, |psi_1> = {format_state(psi_1)}")

    # Step 2: Apply CNOT (control=0, target=1)
    # |psi_2> = CNOT_{1,2} |psi_1>
    psi_2 = apply_gate(psi_1, cnot_on_basis, 0, 1)
    print(f"Step 2: After CNOT(0,1), |psi_2> = {format_state(psi_2)}")

    # Step 3: Apply Toffoli (controls=0,1, target=2)
    # |psi_3> = CCNOT_{1,2,3} |psi_2>
    psi_3 = apply_gate(psi_2, toffoli_on_basis, 0, 1, 2)
    print(f"Step 3: After CCNOT(0,1,2), |psi_3> = {format_state(psi_3)}")

    # Step 4: Apply H to the first qubit again
    # |psi_4> = (H @ I @ I) |psi_3>
    psi_4 = apply_gate(psi_3, hadamard_on_basis, 0)
    print(f"\nStep 4: After H on qubit 0, the final state is:")
    # |psi_4> = 0.5|000> + 0.5|100> + 0.5|011> - 0.5|111>
    final_state_str = "0.5000|000> + 0.5000|011> + 0.5000|100> - 0.5000|111>"
    print(f"|psi_4> = {final_state_str}")
    print("-" * 20)

    # --- Probability Calculation ---
    target_outcome = '100'
    
    # Get the amplitude of the target outcome from the final state
    amplitude = psi_4.get(target_outcome, 0)
    prob = abs(amplitude)**2

    print(f"The amplitude of the outcome |{target_outcome}> in |psi_4> is {amplitude.real:.4f}")
    print(f"The probability of measuring |{target_outcome}> is the squared magnitude of its amplitude.")
    print(f"Probability = |{amplitude.real:.4f}|^2 = {prob:.4f}")
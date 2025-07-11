import math
from collections import defaultdict

def apply_gate(state_dict, gate_func, *qubit_indices):
    """
    Helper function to apply a gate to the state.
    It takes a state dictionary, a gate function, and the qubit indices it acts on.
    """
    new_state = defaultdict(complex)
    for basis_state, amplitude in state_dict.items():
        # The gate function returns a list of (new_basis_state, new_amplitude_factor) tuples
        new_terms = gate_func(basis_state, *qubit_indices)
        for new_basis, amp_factor in new_terms:
            new_state[new_basis] += amplitude * amp_factor
    return dict(new_state)

def hadamard_func(basis_state, qubit_index):
    """Gate logic for Hadamard."""
    sqrt2_inv = 1 / math.sqrt(2)
    state_list = list(basis_state)
    qubit_val = int(state_list[qubit_index])
    
    # |0> -> 1/sqrt(2)*(|0> + |1>)
    state_list[qubit_index] = '0'
    term1_basis = "".join(state_list)
    # |1> -> 1/sqrt(2)*(|0> - |1>)
    state_list[qubit_index] = '1'
    term2_basis = "".join(state_list)
    
    if qubit_val == 0:
        return [(term1_basis, sqrt2_inv), (term2_basis, sqrt2_inv)]
    else: # qubit_val == 1
        return [(term1_basis, sqrt2_inv), (term2_basis, -sqrt2_inv)]

def cnot_func(basis_state, control_index, target_index):
    """Gate logic for CNOT."""
    if basis_state[control_index] == '1':
        state_list = list(basis_state)
        target_val = int(state_list[target_index])
        state_list[target_index] = str(1 - target_val) # Flip target
        return [("".join(state_list), 1.0)]
    else:
        return [(basis_state, 1.0)] # No change

def ccnot_func(basis_state, control1_index, control2_index, target_index):
    """Gate logic for Toffoli (CCNOT)."""
    if basis_state[control1_index] == '1' and basis_state[control2_index] == '1':
        state_list = list(basis_state)
        target_val = int(state_list[target_index])
        state_list[target_index] = str(1 - target_val) # Flip target
        return [("".join(state_list), 1.0)]
    else:
        return [(basis_state, 1.0)] # No change

def main():
    """
    Main function to run the quantum circuit simulation and calculate the probability.
    """
    # Step 0: Initial State |psi_0> = |000>
    psi = {'000': 1.0 + 0.0j}

    # Step 1: Apply Hadamard gate to the first qubit (index 0)
    psi = apply_gate(psi, hadamard_func, 0)

    # Step 2: Apply CNOT gate (control=0, target=1)
    psi = apply_gate(psi, cnot_func, 0, 1)
    
    # Step 3: Apply Toffoli gate (controls=0,1, target=2)
    psi = apply_gate(psi, ccnot_func, 0, 1, 2)
    
    # Step 4: Apply second Hadamard gate to the first qubit
    psi_final = apply_gate(psi, hadamard_func, 0)
    
    # --- Result Analysis ---
    target_state = '100'
    amplitude = psi_final.get(target_state, 0.0j)
    probability = abs(amplitude)**2

    print("The final state |ψ₄> is calculated by applying the quantum gates sequentially.")
    
    # Format the final state string for clean printing
    final_state_terms = []
    # Sorting ensures a consistent output order
    for basis in sorted(psi_final.keys()):
        amp = psi_final[basis]
        # Skip terms with negligible amplitude
        if abs(amp) < 1e-9:
            continue
        
        # Determine the sign for nice printing
        sign = "+ " if amp.real > 0 else "- "
        if len(final_state_terms) == 0:
            sign = "" if amp.real > 0 else "- " # No plus for the first term
            
        final_state_terms.append(f"{sign}{abs(amp.real):.2f}|{basis}>")
        
    print(f"Final State |ψ₄> = {' '.join(final_state_terms)}")

    print(f"\nThe probability of measuring the outcome |{target_state}> is the squared magnitude of its amplitude.")
    # Output the final equation with the calculated numbers
    print(f"P(|{target_state}>) = |(amplitude of |{target_state}>)|^2 = |{amplitude.real:.2f}|^2 = {probability:.2f}")

if __name__ == "__main__":
    main()
    print("\n<<<0.25>>>")
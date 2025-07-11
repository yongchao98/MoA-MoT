import numpy as np

def solve_quantum_circuit():
    """
    Simulates the 3-qubit quantum circuit and calculates the probability of measuring |100>.
    """
    # Step 1: Define the initial state and basic quantum gates

    # Initial state |ψ₀⟩ = |000⟩
    # Basis: |000⟩, |001⟩, |010⟩, |011⟩, |100⟩, |101⟩, |110⟩, |111⟩
    # Index:    0      1      2      3      4      5      6      7
    psi_0 = np.zeros(8, dtype=complex)
    psi_0[0] = 1

    # Single-qubit gates
    I = np.identity(2, dtype=complex)
    H = 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=complex)

    # Step 2: Define the full 3-qubit gate matrices for the circuit

    # Gate for steps 1 & 4: Hadamard on the first qubit
    H1 = np.kron(H, np.kron(I, I))

    # Gate for step 2: CNOT with qubit 1 as control, qubit 2 as target
    # This gate flips the basis states |10x⟩ to |11x⟩ and |11x⟩ to |10x⟩.
    # |100⟩ (idx 4) maps to |110⟩ (idx 6)
    # |101⟩ (idx 5) maps to |111⟩ (idx 7)
    CNOT_1_2 = np.identity(8, dtype=complex)
    CNOT_1_2[[4, 6]] = CNOT_1_2[[6, 4]] # Swap rows 4 and 6
    CNOT_1_2[[5, 7]] = CNOT_1_2[[7, 5]] # Swap rows 5 and 7

    # Gate for step 3: Toffoli (CCNOT) with qubits 1,2 as control, qubit 3 as target
    # This gate flips |110⟩ to |111⟩ and vice versa.
    # |110⟩ (idx 6) maps to |111⟩ (idx 7)
    CCNOT_1_2_3 = np.identity(8, dtype=complex)
    CCNOT_1_2_3[[6, 7]] = CCNOT_1_2_3[[7, 6]] # Swap rows 6 and 7

    # Step 3: Apply the gates in sequence to the state vector

    # After step 1: |ψ₁⟩ = H₁ |ψ₀⟩
    psi_1 = H1 @ psi_0
    
    # After step 2: |ψ₂⟩ = CNOT₁,₂ |ψ₁⟩
    psi_2 = CNOT_1_2 @ psi_1

    # After step 3: |ψ₃⟩ = CCNOT₁,₂,₃ |ψ₂⟩
    psi_3 = CCNOT_1_2_3 @ psi_2

    # After step 4: |ψ₄⟩ = H₁ |ψ₃⟩
    psi_4 = H1 @ psi_3

    # Step 4: Determine the final equation and calculate the probability

    # Build the string for the final state equation
    equation_parts = []
    basis_states = [f'|{i:03b}⟩' for i in range(8)]
    for i, amp in enumerate(psi_4):
        # Check for non-zero amplitudes (within a small tolerance)
        if not np.isclose(amp, 0):
            real_part = np.real(amp)
            sign = "+ " if real_part > 0 else "- "
            
            # Format coefficient to be more readable, e.g., 1/2 instead of 0.5
            val = abs(real_part)
            if np.isclose(val, 0.5):
                coeff_str = "1/2"
            else:
                coeff_str = f"{val:.3f}"
            
            equation_parts.append(f"{sign}{coeff_str}{basis_states[i]}")

    # Format the final equation string
    final_equation_str = " ".join(equation_parts).lstrip('+ ')
    print("The final state of the system is:")
    print(f"|ψ₄⟩ = {final_equation_str}\n")
    
    # The outcome |100⟩ corresponds to index 4 in the state vector
    target_index = 4
    amplitude_100 = psi_4[target_index]
    
    # Probability is the squared magnitude of the amplitude
    probability_100 = np.abs(amplitude_100)**2

    print("The probability of measuring the outcome |100⟩ is calculated as follows:")
    print(f"P(|100⟩) = |⟨100|ψ₄⟩|² = |{np.real(amplitude_100):.3f}|²")
    print(f"P(|100⟩) = {probability_100:.4f}")
    
    return probability_100

# Execute the function and print the final result in the specified format
if __name__ == "__main__":
    result = solve_quantum_circuit()
    # The final answer is 0.25 or 1/4.
    print(f"\n<<<{result:.2f}>>>")
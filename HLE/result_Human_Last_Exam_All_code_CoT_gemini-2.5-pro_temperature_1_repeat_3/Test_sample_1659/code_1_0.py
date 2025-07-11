import numpy as np

def solve_quantum_circuit():
    """
    This function simulates the given 3-qubit quantum circuit and calculates
    the probability of measuring the state |100>.
    """
    # Define single-qubit gates and constants
    I = np.identity(2, dtype=complex)
    H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

    # --- Construct 3-qubit Gate Matrices ---

    # 1. Hadamard on the first qubit (H ⊗ I ⊗ I)
    H1 = np.kron(H, np.kron(I, I))

    # 2. CNOT with qubit 1 as control and 2 as target.
    # This gate flips the second qubit if the first is |1>.
    # It swaps basis states |10z> with |11z>.
    # |100> (4) <-> |110> (6)
    # |101> (5) <-> |111> (7)
    CNOT12 = np.identity(8, dtype=complex)
    CNOT12[[4, 6]] = CNOT12[[6, 4]]
    CNOT12[[5, 7]] = CNOT12[[7, 5]]

    # 3. Toffoli (CCNOT) with qubits 1 and 2 as control, 3 as target.
    # This gate flips the third qubit if the first two are |1>.
    # It swaps basis states |110> with |111>.
    # |110> (6) <-> |111> (7)
    CCNOT123 = np.identity(8, dtype=complex)
    CCNOT123[[6, 7]] = CCNOT123[[7, 6]]

    # --- Simulate the Circuit Step-by-Step ---

    # Initial state: |ψ₀⟩ = |000⟩
    # |000⟩ is the first basis vector in the 8-dim space (index 0).
    psi_0 = np.zeros(8, dtype=complex)
    psi_0[0] = 1

    # Step 1: Apply H to the first qubit
    psi_1 = H1 @ psi_0

    # Step 2: Apply CNOT_1,2
    psi_2 = CNOT12 @ psi_1

    # Step 3: Apply Toffoli gate
    psi_3 = CCNOT123 @ psi_2

    # Step 4: Apply H to the first qubit again
    psi_4 = H1 @ psi_3

    # --- Calculate and Output the Final Probability ---

    # The outcome |100⟩ corresponds to index 4 in the state vector (since 100 in binary is 4).
    index_100 = 4

    # The amplitude is the complex coefficient of the |100⟩ basis state in |ψ₄⟩.
    amplitude_100 = psi_4[index_100]

    # The probability is the squared magnitude of the amplitude.
    probability_100 = np.abs(amplitude_100)**2

    print("--- Quantum Circuit Simulation Result ---")
    print("The final state |ψ₄⟩ has been calculated.")
    print("\nThe probability of measuring the outcome |100⟩ is given by P(|100⟩) = |Amplitude(|100⟩)|².")
    print("\nHere are the numbers for the final equation:")
    print(f"Amplitude(|100⟩) = {amplitude_100.real:.4f}")
    print(f"Probability P(|100⟩) = |{amplitude_100.real:.4f}|^2 = {probability_100:.4f}")

if __name__ == "__main__":
    solve_quantum_circuit()
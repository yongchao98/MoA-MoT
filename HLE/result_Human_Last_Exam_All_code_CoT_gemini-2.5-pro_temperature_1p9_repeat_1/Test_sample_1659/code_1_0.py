import numpy as np

def calculate_quantum_probability():
    """
    This function simulates the given 3-qubit quantum circuit and computes
    the probability of measuring the outcome |100>.
    """

    # --- Define states and gates as matrices ---

    # Single-qubit identity and Hadamard gates
    I = np.identity(2, dtype=complex)
    H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

    # Multi-qubit gate operators for the 3-qubit system
    # H on qubit 1: H ⊗ I ⊗ I
    H1 = np.kron(H, np.kron(I, I))

    # CNOT with control=1, target=2.
    # This is a permutation matrix that swaps |10x> with |11x>.
    # It swaps basis vectors |100> (index 4) ↔ |110> (index 6)
    # and |101> (index 5) ↔ |111> (index 7).
    CNOT12 = np.identity(8, dtype=complex)
    CNOT12[[4, 6]] = CNOT12[[6, 4]]
    CNOT12[[5, 7]] = CNOT12[[7, 5]]

    # CCNOT (Toffoli) with control=1,2, target=3.
    # This is a permutation matrix that swaps |11x>.
    # It swaps basis vectors |110> (index 6) ↔ |111> (index 7).
    CCNOT123 = np.identity(8, dtype=complex)
    CCNOT123[[6, 7]] = CCNOT123[[7, 6]]

    # --- Simulate the quantum circuit ---

    # Initial state |ψ₀⟩ = |000⟩. This is a vector with 1 at index 0.
    psi_0 = np.zeros((8, 1), dtype=complex)
    psi_0[0] = 1

    # Step 1: Apply H to the first qubit
    psi_1 = H1 @ psi_0

    # Step 2: Apply CNOT gate
    psi_2 = CNOT12 @ psi_1

    # Step 3: Apply Toffoli gate
    psi_3 = CCNOT123 @ psi_2

    # Step 4: Apply H to the first qubit again
    psi_4 = H1 @ psi_3

    # --- Calculate Measurement Probability ---

    # The basis state |100> corresponds to index 4 in the state vector
    # (since binary 100 is decimal 4).
    # We extract the complex amplitude for this state.
    amplitude_100 = psi_4[4, 0]

    # The probability is the squared magnitude of the amplitude.
    probability_100 = np.abs(amplitude_100)**2

    # Print the final calculation as requested
    print(f"P(|100⟩) = |{amplitude_100.real:.2f}|^2 = {probability_100:.2f}")

calculate_quantum_probability()
<<<0.25>>>
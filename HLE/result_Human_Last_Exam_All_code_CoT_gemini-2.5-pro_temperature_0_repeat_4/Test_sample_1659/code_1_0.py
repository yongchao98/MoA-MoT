import numpy as np

def solve_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability
    of measuring the state |100>.
    """
    # Define single-qubit gates and identity
    I = np.identity(2, dtype=complex)
    H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

    # Define the initial state |ψ₀⟩ = |000⟩
    # The basis is |000⟩, |001⟩, ..., |111⟩
    psi_0 = np.zeros((8, 1), dtype=complex)
    psi_0[0] = 1

    # Step 1: Apply Hadamard gate to the first qubit: U₁ = H ⊗ I ⊗ I
    U1 = np.kron(H, np.kron(I, I))
    psi_1 = U1 @ psi_0

    # Step 2: Apply CNOT gate with qubit 1 as control and qubit 2 as target.
    # This gate swaps |10z⟩ with |11z⟩.
    # In our basis, this means swapping |100⟩ (4) with |110⟩ (6)
    # and |101⟩ (5) with |111⟩ (7).
    U2 = np.identity(8, dtype=complex)
    # To create the permutation matrix, we swap the corresponding columns of the identity matrix.
    U2[:, [4, 6]] = U2[:, [6, 4]]  # Swaps columns 4 and 6
    U2[:, [5, 7]] = U2[:, [7, 5]]  # Swaps columns 5 and 7
    psi_2 = U2 @ psi_1

    # Step 3: Apply Toffoli gate (CCNOT) with qubits 1,2 as controls and 3 as target.
    # This gate swaps |110⟩ (6) with |111⟩ (7).
    U3 = np.identity(8, dtype=complex)
    U3[:, [6, 7]] = U3[:, [7, 6]]  # Swaps columns 6 and 7
    psi_3 = U3 @ psi_2

    # Step 4: Apply a second Hadamard gate to the first qubit.
    # This is the same operator as in Step 1.
    U4 = U1
    psi_4 = U4 @ psi_3

    # Determine the probability of measuring |100⟩.
    # The basis state |100⟩ corresponds to index 4 (binary 100 = 4).
    amplitude_100 = psi_4[4, 0]
    probability_100 = np.abs(amplitude_100)**2

    # Print the final equation for the probability calculation.
    print("The final state |ψ₄⟩ has been calculated.")
    print("The probability of measuring the outcome |100⟩ is P(|100⟩) = |⟨100|ψ₄⟩|².")
    print("\nThe final equation is:")
    print(f"P(|100⟩) = |{amplitude_100.real:.4f}|^2 = {probability_100:.4f}")

solve_quantum_circuit()
<<<0.25>>>
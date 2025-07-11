import numpy as np

def solve_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability of a specific outcome.
    """
    # --- Step 1: Define single-qubit and multi-qubit gates ---

    # Identity gate
    I = np.identity(2, dtype=complex)
    # Hadamard gate
    H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    # CNOT gate (for 2 qubits)
    CNOT_2q = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 0]], dtype=complex)
    # Toffoli (CCNOT) gate (for 3 qubits)
    CCNOT_3q = np.identity(8, dtype=complex)
    # Swap the |110> and |111> states
    CCNOT_3q[6, 6], CCNOT_3q[7, 7] = 0, 0
    CCNOT_3q[6, 7], CCNOT_3q[7, 6] = 1, 1

    # --- Step 2: Construct the full 8x8 gate matrices ---

    # Hadamard on the first qubit: H ⊗ I ⊗ I
    H1 = np.kron(H, np.kron(I, I))

    # CNOT with q1 as control, q2 as target: CNOT_1,2 ⊗ I
    CNOT12 = np.kron(CNOT_2q, I)

    # Toffoli gate is already defined as an 8x8 matrix
    CCNOT123 = CCNOT_3q

    # --- Step 3: Initialize the state vector ---
    # Initial state |ψ₀⟩ = |000⟩. Index 0 is for |000> (binary 000 = decimal 0).
    psi_0 = np.zeros(8, dtype=complex)
    psi_0[0] = 1

    # --- Step 4: Apply the sequence of gates ---

    # 1. Apply a Hadamard gate to the first qubit
    psi_1 = H1 @ psi_0

    # 2. Apply a CNOT gate with the first qubit as control and the second as target
    psi_2 = CNOT12 @ psi_1

    # 3. Apply a Toffoli gate
    psi_3 = CCNOT123 @ psi_2

    # 4. Apply a second Hadamard gate to the first qubit
    psi_4 = H1 @ psi_3

    # --- Step 5: Calculate the final probability ---

    # The outcome |100⟩ corresponds to index 4 (binary 100 = decimal 4)
    outcome_index = 4
    amplitude = psi_4[outcome_index]
    probability = np.abs(amplitude)**2

    # --- Step 6: Print the result ---
    print("The problem is to find the probability of measuring the state |100>.")
    print("This probability is the squared magnitude of the amplitude of |100> in the final state vector.")
    print("\nCalculation:")
    print(f"Final state amplitude for |100>: {amplitude.real:.4f}")
    print(f"Probability P(|100>) = |{amplitude.real:.4f}|^2 = {probability:.4f}")

# Execute the simulation and print the result
solve_quantum_circuit()
<<<0.25>>>
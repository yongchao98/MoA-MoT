import numpy as np

def solve_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability
    of measuring the state |100>.
    """
    # --- Define initial state and gates ---
    # Initial state |psi_0> = |000>. The basis is ordered as |000>, |001>, ..., |111>.
    # So |000> is the vector [1, 0, 0, 0, 0, 0, 0, 0].
    psi = np.zeros(8, dtype=complex)
    psi[0] = 1

    # Define single-qubit gates
    H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # --- Construct 3-qubit operators for each step ---
    # Step 1 & 4 Operator: Hadamard on qubit 1
    # H_1 = H_q1 ⊗ I_q2 ⊗ I_q3
    H_on_q1 = np.kron(H, np.kron(I, I))

    # Step 2 Operator: CNOT with q1 as control, q2 as target.
    # This gate flips |10x> to |11x> and vice-versa.
    # This corresponds to swapping basis states |100> (index 4) with |110> (index 6)
    # and |101> (index 5) with |111> (index 7).
    CNOT_12 = np.identity(8, dtype=complex)
    CNOT_12[[4, 6], :] = CNOT_12[[6, 4], :]
    CNOT_12[[5, 7], :] = CNOT_12[[7, 5], :]

    # Step 3 Operator: Toffoli (CCNOT) with q1, q2 as controls, q3 as target.
    # This gate flips |11x>, so it swaps |110> with |111>.
    # This corresponds to swapping basis states at index 6 and 7.
    CCNOT_123 = np.identity(8, dtype=complex)
    CCNOT_123[[6, 7], :] = CCNOT_123[[7, 6], :]

    # --- Simulate the circuit by applying gates sequentially ---
    # psi_1 = H_on_q1 @ psi_0
    psi_1 = H_on_q1 @ psi
    # psi_2 = CNOT_12 @ psi_1
    psi_2 = CNOT_12 @ psi_1
    # psi_3 = CCNOT_123 @ psi_2
    psi_3 = CCNOT_123 @ psi_2
    # psi_4 = H_on_q1 @ psi_3
    psi_4 = H_on_q1 @ psi_3

    # --- Determine the probability of measuring |100> ---
    # The outcome |100> corresponds to index 4 (since 100 in binary is 4).
    index_100 = 4
    amplitude_100 = psi_4[index_100]
    probability_100 = np.abs(amplitude_100)**2

    print("The final state is a superposition of |000>, |100>, |011>, and |111>.")
    print("The final state vector is:")
    print(np.round(psi_4, 4))
    
    print("\nThe amplitude for the desired state |100> is the component at index 4.")

    # Print the final equation with its numbers, as requested.
    # The amplitude is 1/2 and the probability is 1/4 or 0.25.
    print("\nFinal Calculation:")
    print(f"P(|100>) = |Amplitude(|100>)|^2 = |{amplitude_100.real:.2f}|^2 = {probability_100:.2f}")

solve_quantum_circuit()
<<<0.25>>>
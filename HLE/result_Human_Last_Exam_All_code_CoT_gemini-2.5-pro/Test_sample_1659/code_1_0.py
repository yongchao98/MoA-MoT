import numpy as np

def solve_quantum_circuit():
    """
    Calculates the probability of measuring the |100> state for a given 3-qubit quantum circuit.
    """
    # --- Define basic quantum objects ---
    # Use complex numbers for all calculations to handle quantum amplitudes
    dtype = np.complex128

    # Qubit basis states
    q0 = np.array([[1], [0]], dtype=dtype)  # |0>
    q1 = np.array([[0], [1]], dtype=dtype)  # |1>

    # Single-qubit gates
    H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=dtype)
    I = np.identity(2, dtype=dtype)

    # Multi-qubit gates
    # CNOT with control on qubit 1, target on qubit 2
    CNOT_12 = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 0]], dtype=dtype)

    # Toffoli (CCNOT) with controls on qubits 1,2 and target on qubit 3
    CCNOT_123 = np.identity(8, dtype=dtype)
    CCNOT_123[6, 6] = 0
    CCNOT_123[7, 7] = 0
    CCNOT_123[6, 7] = 1
    CCNOT_123[7, 6] = 1

    # --- Follow the quantum circuit step-by-step ---

    print("This script calculates the probability of measuring the |100> state.")
    
    # Step 0: Initial state |psi_0> = |000>
    psi_0 = np.kron(q0, np.kron(q0, q0))
    print("Initial State: |psi_0> = |000>\n")

    # Step 1: Apply Hadamard to the first qubit. Operator U1 = H (x) I (x) I
    U1 = np.kron(H, np.kron(I, I))
    psi_1 = np.dot(U1, psi_0)

    # Step 2: Apply CNOT(1,2). Operator U2 = CNOT_12 (x) I
    U2 = np.kron(CNOT_12, I)
    psi_2 = np.dot(U2, psi_1)

    # Step 3: Apply Toffoli CCNOT(1,2,3). Operator U3 = CCNOT_123
    U3 = CCNOT_123
    psi_3 = np.dot(U3, psi_2)

    # Step 4: Apply a second Hadamard to the first qubit. Operator U4 = H (x) I (x) I
    U4 = U1
    psi_4 = np.dot(U4, psi_3)

    # --- Determine the probability of measuring |100> ---
    # The state |100> corresponds to index 4 (binary 100) in the 8-element state vector.
    index_100 = 4

    # The amplitude of the |100> state is the element at this index.
    amplitude_100 = psi_4[index_100][0]

    # The probability is the squared magnitude of the amplitude.
    probability_100 = np.abs(amplitude_100)**2

    # --- Print the final result and the equation ---
    print("The final state |psi_4> is a superposition of four basis states:")
    print("|psi_4> = 0.50|000> + 0.50|100> + 0.50|011> - 0.50|111>")
    print("\nFrom this final state, we can find the probability of measuring |100>.")
    print("\nThe final equation for the probability is:")
    # The instruction: "output each number in the final equation!"
    print(f"P(|100>) = |Amplitude(|100>)|² = |{amplitude_100.real:.2f}|² = {probability_100:.2f}")

if __name__ == '__main__':
    solve_quantum_circuit()
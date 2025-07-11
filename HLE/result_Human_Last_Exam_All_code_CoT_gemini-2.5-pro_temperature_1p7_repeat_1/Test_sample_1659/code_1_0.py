import numpy as np

def solve_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability
    of measuring the state |100>.
    """
    # --- Define Quantum Gates as Matrices ---

    # Single-qubit Identity gate
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    
    # Single-qubit Hadamard gate
    H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

    # 2-qubit CNOT gate (control on the first qubit, target on the second)
    CNOT = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ], dtype=complex)

    # 3-qubit Toffoli (CCNOT) gate (controls on the first two qubits, target on the third)
    CCNOT = np.identity(8, dtype=complex)
    CCNOT[6, 6], CCNOT[7, 7] = 0, 0
    CCNOT[6, 7], CCNOT[7, 6] = 1, 1

    # --- Initial State ---
    # The initial state is |psi_0> = |000>.
    # In vector form, this is [1, 0, 0, 0, 0, 0, 0, 0]^T
    # The basis is ordered |000>, |001>, ..., |111>.
    psi_0 = np.zeros(8, dtype=complex)
    psi_0[0] = 1

    # --- Apply Gates Sequentially ---

    # 1. Apply a Hadamard gate to the first qubit: U = H (x) I (x) I
    U1 = np.kron(H, np.kron(I, I))
    psi_1 = U1 @ psi_0

    # 2. Apply a CNOT gate (control=1st, target=2nd): U = CNOT_1,2 (x) I
    U2 = np.kron(CNOT, I)
    psi_2 = U2 @ psi_1

    # 3. Apply a Toffoli gate (controls=1st,2nd, target=3rd): U = CCNOT_1,2,3
    U3 = CCNOT
    psi_3 = U3 @ psi_2

    # 4. Apply a second Hadamard gate to the first qubit: U = H (x) I (x) I
    U4 = U1
    psi_4 = U4 @ psi_3

    # --- Determine the Probability ---
    # The state |100> corresponds to index 4 in the vector (since 1*2^2 + 0*2^1 + 0*2^0 = 4)
    index_100 = 4
    
    # Get the amplitude of the |100> state from the final state vector.
    amplitude_100 = psi_4[index_100]

    # The probability is the squared magnitude of the amplitude.
    probability_100 = np.abs(amplitude_100)**2

    # --- Print the Final Equation and Result ---
    # The problem requires showing the numbers in the final equation.
    # The equation is P(|100>) = |amplitude|^2 = probability.
    # From calculation, amplitude = 1/2 and probability = 1/4 = 0.25.
    
    print("The final equation for the probability of measuring |100> is:")
    print("P(|100>) = |amplitude|^2")
    print(f"P(|100>) = |{1}/{2}|^2 = {1}/{4} = {probability_100}")

solve_quantum_circuit()
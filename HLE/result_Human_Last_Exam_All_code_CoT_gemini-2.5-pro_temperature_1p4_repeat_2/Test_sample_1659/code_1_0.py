import numpy as np

def run_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the 
    probability of measuring the |100> state.
    """
    
    # --- Define basis states and gates ---
    # Single qubit states
    ket0 = np.array([1, 0], dtype=complex)
    
    # Single qubit gates
    I = np.identity(2, dtype=complex)
    H = (1 / np.sqrt(2)) * np.array([[1, 1], 
                                     [1, -1]], dtype=complex)

    # --- Step 0: Initial State ---
    # The initial state is |psi_0> = |000>
    psi_0 = np.kron(np.kron(ket0, ket0), ket0)
    print("Step 0: Initial State |psi_0> = |000>")
    
    # --- Step 1: Apply H to the first qubit ---
    # Operator is H_1 = H x I x I
    H1_op = np.kron(H, np.kron(I, I))
    psi_1 = H1_op @ psi_0
    print("Step 1: State |psi_1> = (H @ I @ I)|psi_0> = 1/sqrt(2) * (|000> + |100>)")

    # --- Step 2: Apply CNOT(1,2) ---
    # CNOT_1,2 flips the second qubit if the first is 1.
    # It swaps |100> <-> |110> and |101> <-> |111>
    CNOT12_op = np.identity(8, dtype=complex)
    # This matrix swaps basis vector 4 with 6 and 5 with 7
    indices_to_swap = [(4, 6), (5, 7)]
    for i, j in indices_to_swap:
        CNOT12_op[i, i], CNOT12_op[j, j] = 0, 0
        CNOT12_op[i, j], CNOT12_op[j, i] = 1, 1
    psi_2 = CNOT12_op @ psi_1
    print("Step 2: State |psi_2> = CNOT_1,2|psi_1> = 1/sqrt(2) * (|000> + |110>)")
    
    # --- Step 3: Apply Toffoli(1,2,3) ---
    # CCNOT_1,2,3 flips the third qubit if the first two are 1.
    # It swaps |110> <-> |111>
    CCNOT_op = np.identity(8, dtype=complex)
    # This matrix swaps basis vector 6 with 7
    CCNOT_op[6, 6], CCNOT_op[7, 7] = 0, 0
    CCNOT_op[6, 7], CCNOT_op[7, 6] = 1, 1
    psi_3 = CCNOT_op @ psi_2
    print("Step 3: State |psi_3> = CCNOT_1,2,3|psi_2> = 1/sqrt(2) * (|000> + |111>)")
    
    # --- Step 4: Apply H to the first qubit again ---
    psi_4 = H1_op @ psi_3
    print("Step 4: Final State |psi_4> = (H @ I @ I)|psi_3> = 1/2 * (|000> + |100> + |011> - |111>)")
    print("-" * 20)

    # --- Determine the probability of measuring |100> ---
    # The state |100> corresponds to index 4 in the state vector (binary 100 = 4)
    index_100 = 4
    amplitude_100 = psi_4[index_100]
    probability_100 = np.abs(amplitude_100)**2

    print("The final state can be written as:")
    print("|psi_4> = (1/2)|000> + (1/2)|011> + (1/2)|100> - (1/2)|111>")
    
    print("\nThe amplitude for the outcome |100> is its coefficient in the final state.")
    print(f"Amplitude( |100> ) = {amplitude_100.real:.2f}")
    
    print("\nThe probability is the squared magnitude of this amplitude.")
    print(f"P(|100>) = |{amplitude_100.real:.2f}|^2 = {probability_100:.2f}")

    print("\nFinal calculation expressed with fractions:")
    # Printing each number in the final equation as requested
    print("P(|100>) = |(1/2)|^2 = 1/4 = 0.25")

if __name__ == '__main__':
    run_quantum_circuit()
<<<0.25>>>
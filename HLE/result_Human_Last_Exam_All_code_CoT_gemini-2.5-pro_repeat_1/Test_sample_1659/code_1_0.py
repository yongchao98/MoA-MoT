import numpy as np

def run_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability of measuring |100>.
    """
    # Define single-qubit gates
    I = np.eye(2, dtype=complex)
    H = (1 / np.sqrt(2)) * np.array([[1, 1], 
                                     [1, -1]], dtype=complex)

    # --- Step 1: Initial State and First Gate ---
    # Initial state |psi_0> = |000>
    # The basis is |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
    psi_0 = np.zeros(8, dtype=complex)
    psi_0[0] = 1

    # Apply H to the first qubit: H_1 = H x I x I
    H_1 = np.kron(H, np.kron(I, I))
    psi_1 = H_1 @ psi_0

    # --- Step 2: CNOT Gate ---
    # Apply CNOT with qubit 1 as control, qubit 2 as target.
    # This flips the 2nd qubit if the 1st is 1.
    # It swaps |10x> with |11x>.
    CNOT_12 = np.eye(8, dtype=complex)
    # Swap |100> (4) and |110> (6)
    CNOT_12[[4, 6], :] = CNOT_12[[6, 4], :]
    # Swap |101> (5) and |111> (7)
    CNOT_12[[5, 7], :] = CNOT_12[[7, 5], :]
    psi_2 = CNOT_12 @ psi_1

    # --- Step 3: Toffoli Gate ---
    # Apply Toffoli (CCNOT) with qubits 1,2 as controls, qubit 3 as target.
    # This flips the 3rd qubit if the first two are 1.
    # It swaps |110> with |111>.
    CCNOT_123 = np.eye(8, dtype=complex)
    # Swap |110> (6) and |111> (7)
    CCNOT_123[[6, 7], :] = CCNOT_123[[7, 6], :]
    psi_3 = CCNOT_123 @ psi_2
    
    # --- Step 4: Second Hadamard Gate ---
    # Apply H to the first qubit again.
    psi_4 = H_1 @ psi_3

    # --- Determine the probability of measuring |100> ---
    # The state |100> corresponds to index 4 in the state vector.
    index_100 = 4
    amplitude_100 = psi_4[index_100]
    
    # Probability is the squared magnitude of the amplitude.
    probability_100 = np.abs(amplitude_100)**2

    # Print the final calculation as requested
    # The format shows each number in the final equation.
    print(f"The final state vector is:\n{np.round(psi_4, 3)}\n")
    print("The probability of measuring the state |100> is the squared magnitude of its amplitude.")
    # Use .real because the amplitude in this case is a real number.
    print(f"P(|100>) = |{amplitude_100.real:.4f}|^2 = {probability_100:.4f}")

if __name__ == '__main__':
    run_quantum_circuit()
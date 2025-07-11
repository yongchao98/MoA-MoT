import numpy as np

def ket(state_str):
    """Creates a 3-qubit basis state vector from a binary string like '100'."""
    idx = int(state_str, 2)
    vec = np.zeros(8, dtype=complex)
    vec[idx] = 1
    return vec

def main():
    """
    Simulates the 3-qubit quantum circuit and calculates the probability of measuring |100>.
    """
    # --- Define Single-Qubit Gates ---
    H = (1 / np.sqrt(2)) * np.array([[1, 1],
                                     [1, -1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # --- Define 3-Qubit Circuit Operators ---
    # Hadamard on qubit 1
    H1 = np.kron(H, np.kron(I, I))

    # CNOT with qubit 1 as control and qubit 2 as target
    CNOT12 = np.zeros((8, 8), dtype=complex)
    for i in range(8):
        q1 = (i >> 2) & 1 # Control bit
        q2 = (i >> 1) & 1 # Target bit
        q3 = (i >> 0) & 1
        
        target_idx = i
        if q1 == 1:
            # If control is 1, flip target q2
            flipped_q2 = 1 - q2
            target_idx = (q1 << 2) | (flipped_q2 << 1) | q3
        
        CNOT12[target_idx, i] = 1
        
    # Toffoli (CCNOT) with qubits 1, 2 as control and 3 as target
    # This gate flips |110> to |111> and vice-versa.
    CCNOT = np.identity(8, dtype=complex)
    # The transformation CCNOT |a,b,c> = |a,b, c XOR ab>
    # means it only affects |110> and |111>
    # CCNOT |110> -> |111> and CCNOT |111> -> |110>
    # which is a swap of the corresponding basis vectors
    idx_110 = int('110', 2)
    idx_111 = int('111', 2)
    CCNOT[:, [idx_110, idx_111]] = CCNOT[:, [idx_111, idx_110]]


    # --- Run the Simulation ---
    # 1. Initial state |psi_0> = |000>
    psi_0 = ket('000')

    # 2. Apply gates in sequence
    # psi_1 = H_1 |psi_0>
    # psi_2 = CNOT_{1,2} |psi_1>
    # psi_3 = CCNOT_{1,2,3} |psi_2>
    # psi_4 = H_1 |psi_3>
    
    # We can combine the matrix multiplications
    U_circuit = H1 @ CCNOT @ CNOT12 @ H1
    psi_final = U_circuit @ psi_0

    # --- Determine the Probability ---
    # Find the amplitude of the |100> state.
    # The basis is |q1,q2,q3>, so |100> corresponds to decimal 4.
    index_100 = int('100', 2)
    amplitude_100 = psi_final[index_100]

    # The probability is the squared magnitude of the amplitude.
    probability_100 = np.abs(amplitude_100)**2

    # --- Print the Result ---
    # The user requested the numbers in the final equation to be printed.
    print("The final calculation for the probability of measuring |100> is:")
    print(f"P(|100>) = |amplitude|^2 = |{amplitude_100.real:.4f}|^2 = {probability_100:.4f}")

if __name__ == "__main__":
    main()
import numpy as np

def state_to_string(vec, n_qubits=3, tol=1e-9):
    """Converts a state vector to its string representation in ket notation."""
    terms = []
    # Find the first non-zero term to handle the initial sign
    first_term_found = False
    for i, amp in enumerate(vec):
        if abs(amp) > tol:
            # Get the basis state in ket notation, e.g., |001>
            ket = f"|{i:0{n_qubits}b}>"

            # Determine the sign for the term
            is_negative = amp.real < -tol or (abs(amp.real) < tol and amp.imag < -tol)
            sign = " - " if is_negative else " + "
            
            # The first term does not need a leading "+" sign
            if not first_term_found:
                sign = "-" if is_negative else ""
                first_term_found = True

            # Get the absolute value of the amplitude
            val = np.abs(amp)
            
            # Format the amplitude value
            amp_str = ""
            if not np.isclose(val, 1.0):
                amp_str = f"{val:.3f}"

            terms.append(f"{sign}{amp_str}{ket}")
            
    if not terms:
        return "0"
    return "".join(terms).strip()

def run_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability of measuring |100>.
    """
    n_qubits = 3
    dim = 2**n_qubits

    # Step 0: Initial State |psi_0> = |000>
    # The basis states are ordered from |000> (index 0) to |111> (index 7)
    psi_0 = np.zeros(dim, dtype=complex)
    psi_0[0] = 1.0
    print("Initial State |psi_0> =", state_to_string(psi_0, n_qubits))
    print("-" * 30)

    # --- Define Quantum Gates ---
    # Single-qubit gates
    H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # H gate on the first qubit
    H_1 = np.kron(H, np.kron(I, I))

    # CNOT gate with q1 as control, q2 as target
    # This swaps |10x> with |11x>
    # |100> (idx 4) <-> |110> (idx 6)
    # |101> (idx 5) <-> |111> (idx 7)
    CNOT_1_2 = np.identity(dim, dtype=complex)
    CNOT_1_2[[4, 6], :] = CNOT_1_2[[6, 4], :]
    CNOT_1_2[[5, 7], :] = CNOT_1_2[[7, 5], :]

    # Toffoli (CCNOT) gate with q1,q2 as controls, q3 as target
    # This swaps |11x> with |11(x XOR 1)>
    # |110> (idx 6) <-> |111> (idx 7)
    CCNOT_1_2_3 = np.identity(dim, dtype=complex)
    CCNOT_1_2_3[[6, 7], :] = CCNOT_1_2_3[[7, 6], :]

    # --- Simulate the Circuit ---
    # Step 1: Apply Hadamard gate to the first qubit
    psi_1 = H_1 @ psi_0
    print("Step 1: Apply H to Qubit 1")
    print("|psi_1> = (H*I*I)|psi_0> =", state_to_string(psi_1, n_qubits))
    print("-" * 30)

    # Step 2: Apply CNOT gate
    psi_2 = CNOT_1_2 @ psi_1
    print("Step 2: Apply CNOT(1,2)")
    print("|psi_2> = CNOT_1,2 |psi_1> =", state_to_string(psi_2, n_qubits))
    print("-" * 30)

    # Step 3: Apply Toffoli gate
    psi_3 = CCNOT_1_2_3 @ psi_2
    print("Step 3: Apply CCNOT(1,2,3)")
    print("|psi_3> = CCNOT_1,2,3 |psi_2> =", state_to_string(psi_3, n_qubits))
    print("-" * 30)

    # Step 4: Apply a second Hadamard gate to the first qubit
    psi_4 = H_1 @ psi_3
    print("Step 4: Apply H to Qubit 1 again")
    print("Final state |psi_4> =", state_to_string(psi_4, n_qubits))
    print("=" * 30)
    
    # --- Calculate the Probability ---
    # The state |100> corresponds to index 4 (since binary 100 is decimal 4)
    index_100 = 4
    amplitude = psi_4[index_100]
    probability = np.abs(amplitude)**2

    print("\nCalculating the probability of measuring state |100>:")
    print("The probability is the squared magnitude of the amplitude of |100>.")
    print(f"The final state is: |psi_4> = {state_to_string(psi_4, n_qubits)}")
    print(f"The amplitude for the component |100> is {amplitude.real:.3f}.")
    print("\nThe final probability calculation is:")
    print(f"P(|100>) = |Amplitude(|100>)|^2 = |{amplitude.real:.3f}|^2 = {probability:.3f}")

if __name__ == '__main__':
    run_quantum_circuit()
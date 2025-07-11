import numpy as np

def solve_quantum_circuit():
    """
    Solves the given 3-qubit quantum circuit problem.
    """
    # Define single qubit states and gates
    q0 = np.array([[1], [0]], dtype=complex)  # |0>
    q1 = np.array([[0], [1]], dtype=complex)  # |1>
    I = np.identity(2, dtype=complex)
    H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

    # Step 0: Initial State |psi_0> = |000>
    # Basis order is |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
    psi_0 = np.kron(np.kron(q0, q0), q0)

    # Step 1: Apply Hadamard gate to the first qubit: |psi_1> = (H x I x I)|psi_0>
    G1 = np.kron(np.kron(H, I), I)
    psi_1 = np.dot(G1, psi_0)

    # Step 2: Apply CNOT gate (q1 as control, q2 as target): |psi_2> = CNOT_12 |psi_1>
    # This gate flips the second qubit if the first is |1>.
    # It swaps |10x> with |11x>. In our basis, this means swapping
    # |100> (idx 4) with |110> (idx 6) and |101> (idx 5) with |111> (idx 7).
    CNOT_12 = np.identity(8, dtype=complex)
    CNOT_12[4, 4], CNOT_12[6, 6] = 0, 0
    CNOT_12[4, 6], CNOT_12[6, 4] = 1, 1
    CNOT_12[5, 5], CNOT_12[7, 7] = 0, 0
    CNOT_12[5, 7], CNOT_12[7, 5] = 1, 1
    psi_2 = np.dot(CNOT_12, psi_1)

    # Step 3: Apply Toffoli gate (q1,q2 as control, q3 as target): |psi_3> = CCNOT_123 |psi_2>
    # This gate flips the third qubit if the first two are |11>.
    # It swaps |110> (idx 6) with |111> (idx 7).
    CCNOT_123 = np.identity(8, dtype=complex)
    CCNOT_123[6, 6], CCNOT_123[7, 7] = 0, 0
    CCNOT_123[6, 7], CCNOT_123[7, 6] = 1, 1
    psi_3 = np.dot(CCNOT_123, psi_2)

    # Step 4: Apply a second Hadamard gate to the first qubit: |psi_4> = (H x I x I)|psi_3>
    G4 = G1
    psi_4 = np.dot(G4, psi_3)

    # Find the probability of measuring |100>
    # The basis state |100> corresponds to index 4.
    amplitude_100 = psi_4[4, 0]
    probability_100 = np.abs(amplitude_100)**2

    print("The final state |psi_4> is:")
    final_state_parts = []
    basis_states = [f"|{i:03b}>" for i in range(8)]
    for i, amp in enumerate(psi_4.flatten()):
        if not np.isclose(amp, 0):
            # Format complex numbers nicely
            if np.isclose(amp.imag, 0):
                amp_str = f"{amp.real:.2f}"
            else:
                amp_str = f"{amp:.2f}"
            final_state_parts.append(f"({amp_str}) {basis_states[i]}")
    print(" + ".join(final_state_parts).replace("+ (-", "- ("))

    print("\nTo find the probability of measuring the outcome |100>:")
    print("1. Find the amplitude of the |100> component in the final state.")
    print(f"   The amplitude is c_100 = {amplitude_100.real:.2f}")

    print("2. The probability is the squared magnitude of this amplitude.")
    print(f"   P(|100>) = |c_100|^2 = |{amplitude_100.real:.2f}|^2 = {probability_100:.2f}")

solve_quantum_circuit()
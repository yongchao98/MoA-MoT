import numpy as np

def solve_quantum_circuit():
    """
    Calculates the final state of a 3-qubit system and the probability of a specific outcome.
    """
    # --- 1. Define single-qubit states and gates ---
    ket0 = np.array([1, 0], dtype=complex)
    ket1 = np.array([0, 1], dtype=complex)
    
    # Hadamard Gate
    H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
    # Identity Gate
    I = np.identity(2, dtype=complex)

    # --- 2. Define the initial 3-qubit state ---
    # |psi_0> = |000> = |0> kronecker |0> kronecker |0>
    psi_0 = np.kron(ket0, np.kron(ket0, ket0))

    # --- 3. Construct the 3-qubit gate matrices ---
    # Hadamard on the first qubit
    H1 = np.kron(H, np.kron(I, I))

    # CNOT with qubit 1 as control, qubit 2 as target
    # This gate flips the state |10x> to |11x> and vice-versa.
    # It corresponds to swapping basis states |100> (index 4) with |110> (index 6)
    # and |101> (index 5) with |111> (index 7).
    CNOT12 = np.identity(8, dtype=complex)
    CNOT12[[4, 6], :] = CNOT12[[6, 4], :] # Swap rows for |100> and |110>
    CNOT12[[5, 7], :] = CNOT12[[7, 5], :] # Swap rows for |101> and |111>

    # Toffoli (CCNOT) with qubits 1,2 as control, qubit 3 as target
    # This gate flips the state |11c> to |11(1-c)>.
    # It corresponds to swapping basis states |110> (index 6) with |111> (index 7).
    CCNOT123 = np.identity(8, dtype=complex)
    CCNOT123[[6, 7], :] = CCNOT123[[7, 6], :] # Swap rows for |110> and |111>

    # --- 4. Simulate the quantum circuit ---
    # Step 1: Apply Hadamard to the first qubit
    psi_1 = H1 @ psi_0
    
    # Step 2: Apply CNOT gate
    psi_2 = CNOT12 @ psi_1
    
    # Step 3: Apply Toffoli gate
    psi_3 = CCNOT123 @ psi_2
    
    # Step 4: Apply second Hadamard gate
    psi_4 = H1 @ psi_3

    # --- 5. Analyze the final state and calculate probability ---
    print("The final state |psi_4> is a superposition of four basis states.")
    
    # Build and print the equation for the final state
    basis_kets = [f"|{i:03b}>" for i in range(8)]
    terms = []
    for i, amp in enumerate(psi_4):
        if not np.isclose(np.abs(amp), 0):
            # The amplitudes are real in this specific problem
            terms.append(f"({amp.real:.2f}){basis_kets[i]}")
    
    # Join terms with '+' and clean up for negative signs
    final_equation_str = " + ".join(terms).replace("+ (-", "- (")
    print(f"Final equation: |psi_4> = {final_equation_str}\n")

    # The basis state |100> corresponds to index 4 (since 100 in binary is 4 in decimal)
    index_100 = 4
    amplitude_100 = psi_4[index_100]
    
    # Probability is the squared magnitude of the amplitude
    probability_100 = np.abs(amplitude_100)**2

    print(f"The amplitude of measuring the state |100> is {amplitude_100.real:.2f}.")
    print(f"The probability is the squared magnitude of the amplitude:")
    print(f"P(|100>) = |{amplitude_100.real:.2f}|^2 = {probability_100:.2f}")

# Execute the simulation
solve_quantum_circuit()
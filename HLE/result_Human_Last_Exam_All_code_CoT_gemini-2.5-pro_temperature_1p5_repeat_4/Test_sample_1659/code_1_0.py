import numpy as np

def run_quantum_circuit():
    """
    Simulates the given 3-qubit quantum circuit and calculates the probability of measuring |100>.
    """
    # Number of qubits
    n_qubits = 3
    # Total number of states in the system
    num_states = 2**n_qubits

    # --- 1. Define Initial State ---
    # The initial state is |000>. The index for ket |abc> is 4a + 2b + c.
    # So, |000> is at index 0.
    psi_0 = np.zeros(num_states, dtype=complex)
    psi_0[0] = 1.0
    print("Initial state |psi_0> = |000>")
    
    # --- 2. Define Quantum Gates as Matrices ---
    
    # Single-qubit gates
    H = (1 / np.sqrt(2)) * np.array([[1, 1], 
                                     [1, -1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # H gate on the first qubit (H x I x I)
    H1 = np.kron(H, np.kron(I, I))

    # CNOT gate with qubit 1 as control, qubit 2 as target
    # This swaps |10c> with |11c>
    # |100> (idx 4) <-> |110> (idx 6)
    # |101> (idx 5) <-> |111> (idx 7)
    CNOT12 = np.identity(num_states, dtype=complex)
    CNOT12[[4, 6]] = CNOT12[[6, 4]]
    CNOT12[[5, 7]] = CNOT12[[7, 5]]

    # Toffoli (CCNOT) gate with qubits 1, 2 as control, qubit 3 as target
    # This swaps |110> with |111>
    # |110> (idx 6) <-> |111> (idx 7)
    CCNOT123 = np.identity(num_states, dtype=complex)
    CCNOT123[[6, 7]] = CCNOT123[[7, 6]]

    # --- 3. Apply Gates Sequentially ---

    # Step 1: Apply Hadamard to the first qubit
    psi_1 = H1 @ psi_0

    # Step 2: Apply CNOT(1,2)
    psi_2 = CNOT12 @ psi_1
    
    # Step 3: Apply Toffoli(1,2,3)
    psi_3 = CCNOT123 @ psi_2
    
    # Step 4: Apply second Hadamard to the first qubit
    psi_4 = H1 @ psi_3

    # --- 4. Determine Final State and Probability ---

    print("\nFinal state |psi_4> calculation:")
    
    equation_parts = []
    for i in range(num_states):
        amplitude = psi_4[i]
        if not np.isclose(amplitude, 0):
            ket = f"|{i:03b}>"
            if not equation_parts: # First term
                 # Handle negative sign for the very first term
                if np.isclose(amplitude.real, -1.0):
                    sign = "-"
                    amp_str = ""
                else:
                    sign = "" if amplitude.real > 0 else "-"
                    amp_str = f"{abs(amplitude.real):.4f}"
                equation_parts.append(f"{sign}{amp_str}{ket}")
            else: # Subsequent terms
                sign = " + " if amplitude.real > 0 else " - "
                amp_str = f"{abs(amplitude.real):.4f}"
                equation_parts.append(f"{sign}{amp_str}{ket}")

    print(f"|psi_4> = {''.join(equation_parts)}")
    
    # The state to measure is |100>
    target_state_ket = '100'
    target_index = int(target_state_ket, 2)
    
    # Get the amplitude of the target state
    amplitude = psi_4[target_index]
    
    # Probability is the squared magnitude of the amplitude
    probability = np.abs(amplitude)**2
    
    print(f"\nThe probability of measuring the state |{target_state_ket}> is given by the squared magnitude of its amplitude.")
    print(f"Final equation for the probability:")
    print(f"P(|{target_state_ket}>) = |Amplitude({target_state_ket})|^2 = |{amplitude.real:.4f}|^2 = {probability:.4f}")

    return probability

if __name__ == '__main__':
    final_probability = run_quantum_circuit()
    # The final required output format
    # print(f"\n<<<{final_probability:.2f}>>>")
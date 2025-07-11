import numpy as np

def calculate_quantum_probability():
    """
    Calculates the probability of measuring the |100> state in a 3-qubit system
    after a sequence of specified quantum gate operations.
    """
    # --- 1. Define single-qubit elements ---
    # States |0> and |1> as column vectors
    q0 = np.array([[1], [0]], dtype=complex)
    
    # Hadamard gate (H) and Identity gate (I)
    H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # --- 2. Setup the 3-qubit system ---
    # Initial state |ψ₀⟩ = |000⟩ = |0⟩ ⊗ |0⟩ ⊗ |0⟩
    # We use np.kron for the tensor product.
    psi_0 = np.kron(q0, np.kron(q0, q0))

    # --- 3. Define the 3-qubit gate operators (as 8x8 matrices) ---
    # Operator for Step 1 & 4: Hadamard on qubit 1 (H ⊗ I ⊗ I)
    H_on_q1 = np.kron(H, np.kron(I, I))

    # Operator for Step 2: CNOT with control=1, target=2.
    # This gate flips the target qubit if the control is 1.
    # It swaps the basis states |100> (index 4) with |110> (index 6)
    # and |101> (index 5) with |111> (index 7).
    CNOT_12 = np.identity(8, dtype=complex)
    CNOT_12[[4, 6],:] = CNOT_12[[6, 4],:]
    CNOT_12[[5, 7],:] = CNOT_12[[7, 5],:]
    
    # Operator for Step 3: Toffoli (CCNOT) with controls=1,2, target=3.
    # This gate flips the target if both controls are 1.
    # It swaps the basis states |110> (index 6) and |111> (index 7).
    CCNOT_123 = np.identity(8, dtype=complex)
    CCNOT_123[[6, 7],:] = CCNOT_123[[7, 6],:]
    
    # --- 4. Apply the quantum circuit gates sequentially ---
    # Use matrix multiplication (@) to apply each gate to the state vector.
    psi_1 = H_on_q1 @ psi_0
    psi_2 = CNOT_12 @ psi_1
    psi_3 = CCNOT_123 @ psi_2
    psi_4_final = H_on_q1 @ psi_3

    # --- 5. Determine the probability of measuring |100> ---
    # The basis state |100> (binary) corresponds to index 4 (decimal) in the state vector.
    index_100 = 4
    amplitude_100 = psi_4_final[index_100, 0] # Extract the scalar amplitude
    
    # Probability is the squared magnitude of the amplitude.
    prob_100 = np.abs(amplitude_100)**2
    
    # --- 6. Print the final calculation clearly ---
    final_amp_real = amplitude_100.real
    
    print("This script calculates the probability of measuring the state |100>.")
    print("--- Final Calculation ---")
    print("The amplitude 'c' of the |100> state in the final state vector is found.")
    print(f"c = {final_amp_real:.2f}")
    
    print("\nThe probability P is the squared magnitude of this amplitude: P = |c|²")
    print(f"P(|100>) = |{final_amp_real:.2f}|^2 = {prob_100:.2f}")

# Execute the function
calculate_quantum_probability()
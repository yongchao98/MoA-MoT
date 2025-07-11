import numpy as np

# 1. Initialize States and Gates

# Single-qubit gates
I = np.array([[1, 0], [0, 1]], dtype=complex)
H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

# 2. Construct 3-qubit gate matrices

# Hadamard on the first qubit (H ⊗ I ⊗ I)
H1 = np.kron(H, np.kron(I, I))

# CNOT with qubit 1 as control and 2 as target.
# It maps basis |10x> to |11x> and vice-versa.
# This corresponds to swapping basis vectors |100>(4) with |110>(6)
# and |101>(5) with |111>(7).
CNOT_12 = np.identity(8, dtype=complex)
CNOT_12[[4, 6]] = CNOT_12[[6, 4]] # Swap rows for |100> and |110>
CNOT_12[[5, 7]] = CNOT_12[[7, 5]] # Swap rows for |101> and |111>

# Toffoli (CCNOT) with qubits 1,2 as control and 3 as target.
# It maps basis |110> to |111> and vice-versa.
# This corresponds to swapping basis vectors |110>(6) and |111>(7).
CCNOT_123 = np.identity(8, dtype=complex)
CCNOT_123[[6, 7]] = CCNOT_123[[7, 6]] # Swap rows for |110> and |111>


# 3. Simulate the Quantum Circuit

# Initial state |ψ₀⟩ = |000⟩
# The basis is ordered |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
psi_0 = np.zeros(8, dtype=complex)
psi_0[0] = 1

# Step 1: Apply H to the first qubit
psi_1 = H1 @ psi_0

# Step 2: Apply CNOT(1,2)
psi_2 = CNOT_12 @ psi_1

# Step 3: Apply CCNOT(1,2,3)
psi_3 = CCNOT_123 @ psi_2

# Step 4: Apply H to the first qubit again
psi_4 = H1 @ psi_3

# 4. Determine Probability

# The basis state |100> corresponds to index 4 in the state vector.
index_100 = 4
amplitude_100 = psi_4[index_100]

# Probability is the squared magnitude of the amplitude.
probability_100 = np.abs(amplitude_100)**2

# 5. Print the Final Result
print("The final state vector |ψ₄⟩ is:")
print(np.round(psi_4, 4))
print("\nThe probability of measuring the state |100> is P(|100>) = |<100|ψ₄>|².")
print(f"The amplitude for |100> is <100|ψ₄> = {amplitude_100.real:.4f} + {amplitude_100.imag:.4f}j")
print("\nFinal calculation:")
# The format required: output each number in the final equation.
print(f"P(|100>) = |{amplitude_100:.4f}|² = {probability_100:.4f}")

# Sanity check: verify total probability is 1
total_prob = np.sum(np.abs(psi_4)**2)
print(f"\nVerification: Total probability of all states is {total_prob:.2f}")

<<<0.25>>>
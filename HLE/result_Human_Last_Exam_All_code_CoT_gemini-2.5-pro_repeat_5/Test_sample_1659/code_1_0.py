import numpy as np

# This script simulates a 3-qubit quantum circuit to find the probability of a specific outcome.

# --- Define single-qubit gates and states ---
# Basis states |0> and |1>
q0 = np.array([1, 0], dtype=complex)
q1 = np.array([0, 1], dtype=complex)

# Hadamard gate
H = (1 / np.sqrt(2)) * np.array([[1, 1],
                                 [1, -1]], dtype=complex)
# Identity gate
I = np.identity(2, dtype=complex)

# --- Define 3-qubit gate matrices ---
# H on qubit 1: H ⊗ I ⊗ I
H1 = np.kron(H, np.kron(I, I))

# CNOT with control on qubit 1, target on qubit 2
# This swaps the amplitudes of |10x> and |11x> states.
# Basis order: |000>,|001>,|010>,|011>,|100>,|101>,|110>,|111>
# Indices:      0,    1,    2,    3,    4,    5,    6,    7
CNOT_12 = np.identity(8, dtype=complex)
# Swap |100> (idx 4) with |110> (idx 6)
CNOT_12[[4, 6]] = CNOT_12[[6, 4]]
# Swap |101> (idx 5) with |111> (idx 7)
CNOT_12[[5, 7]] = CNOT_12[[7, 5]]

# Toffoli (CCNOT) with controls on qubits 1 & 2, target on qubit 3
# This swaps the amplitudes of |110> and |111> states.
CCNOT_123 = np.identity(8, dtype=complex)
# Swap |110> (idx 6) with |111> (idx 7)
CCNOT_123[[6, 7]] = CCNOT_123[[7, 6]]

# --- Simulate the quantum circuit ---

# Initial state: |ψ₀⟩ = |000⟩
psi_0 = np.kron(q0, np.kron(q0, q0))

# Step 1: Apply H to the first qubit
# |ψ₁⟩ = (H ⊗ I ⊗ I) |ψ₀⟩
psi_1 = H1 @ psi_0

# Step 2: Apply CNOT with control on qubit 1, target on qubit 2
# |ψ₂⟩ = CNOT₁,₂ |ψ₁⟩
psi_2 = CNOT_12 @ psi_1

# Step 3: Apply Toffoli gate
# |ψ₃⟩ = CCNOT₁,₂,₃ |ψ₂⟩
psi_3 = CCNOT_123 @ psi_2

# Step 4: Apply H to the first qubit again
# |ψ₄⟩ = (H ⊗ I ⊗ I) |ψ₃⟩
psi_4 = H1 @ psi_3

# --- Determine the probability of measuring |100> ---

# The basis state |100> corresponds to index 4 in the state vector (binary 100 = 4)
index_100 = 4
amplitude_100 = psi_4[index_100]

# Probability is the squared magnitude of the amplitude
probability_100 = np.abs(amplitude_100)**2

# --- Print the final calculation ---
print("The final state |ψ₄⟩ results from the sequence of gate operations.")
print("We want to find the probability of measuring the outcome |100⟩.")
print("This probability is the squared magnitude of the amplitude of |100⟩ in the final state.")
print("\nFinal Calculation:")
print(f"Amplitude of |100⟩: c_100 = {amplitude_100.real:.4f}")
print(f"Probability P(|100⟩) = |c_100|² = |{amplitude_100.real:.4f}|² = {probability_100:.4f}")

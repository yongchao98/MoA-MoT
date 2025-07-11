import numpy as np

# --- 1. Define single-qubit gates and initial state ---

# Hadamard gate
H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=complex)
# Identity gate
I = np.eye(2, dtype=complex)

# Initial state |ψ₀⟩ = |000⟩
# The basis is ordered as |000⟩, |001⟩, |010⟩, |011⟩, |100⟩, |101⟩, |110⟩, |111⟩
# Index for |q₁q₂q₃⟩ is 4*q₁ + 2*q₂ + 1*q₃
psi_0 = np.zeros(8, dtype=complex)
psi_0[0] = 1
print(f"Initial state |ψ₀⟩ = |000⟩")

# --- 2. Define 3-qubit circuit operations ---

# Step 1: Hadamard on the first qubit
# Operator is H ⊗ I ⊗ I
H1 = np.kron(H, np.kron(I, I))

# Step 2: CNOT with qubit 1 as control, qubit 2 as target
# This gate flips the second qubit if the first is |1⟩.
# It swaps |10x⟩ with |11x⟩.
# |100⟩ (idx 4) ↔ |110⟩ (idx 6)
# |101⟩ (idx 5) ↔ |111⟩ (idx 7)
CNOT12 = np.eye(8, dtype=complex)
CNOT12[[4, 6], :] = CNOT12[[6, 4], :]
CNOT12[[5, 7], :] = CNOT12[[7, 5], :]

# Step 3: Toffoli (CCNOT) with qubits 1,2 as controls, qubit 3 as target
# This gate flips the third qubit if the first two are |1⟩.
# It swaps |110⟩ with |111⟩.
# |110⟩ (idx 6) ↔ |111⟩ (idx 7)
CCNOT123 = np.eye(8, dtype=complex)
CCNOT123[[6, 7], :] = CCNOT123[[7, 6], :]

# --- 3. Simulate the quantum circuit ---

# Apply first Hadamard gate
psi_1 = H1 @ psi_0

# Apply CNOT gate
psi_2 = CNOT12 @ psi_1

# Apply Toffoli gate
psi_3 = CCNOT123 @ psi_2

# Apply second Hadamard gate
psi_4 = H1 @ psi_3

# --- 4. Calculate and print the result ---

# The target outcome is |100⟩, which corresponds to index 4
target_index = 4
amplitude = psi_4[target_index]
probability = np.abs(amplitude)**2

# Print the final equation for the probability
print(f"\nThe final state |ψ₄⟩ has been calculated.")
print(f"The amplitude of the target state |100⟩ is {amplitude.real:.4f}")
print("\nThe probability of measuring the state |100⟩ is given by the squared magnitude of its amplitude.")
# The prompt requires outputting each number in the final equation.
print(f"P(|100⟩) = |c_100|² = |{amplitude.real:.2f}|² = {probability:.4f}")

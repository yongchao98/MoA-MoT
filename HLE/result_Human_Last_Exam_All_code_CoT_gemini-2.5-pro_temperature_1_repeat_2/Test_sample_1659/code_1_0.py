import numpy as np

# --- Define single-qubit gates and states ---
I = np.eye(2, dtype=complex)
H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

# --- Define 3-qubit operators ---
# H on qubit 1: H_1 = H_q1 x I_q2 x I_q3
H1 = np.kron(H, np.kron(I, I))

# CNOT with control=1, target=2
# This can be constructed by finding its permutation matrix
# It swaps |10x> with |11x>
CNOT12 = np.eye(8, dtype=complex)
# Swap |100> (4) and |110> (6)
CNOT12[[4, 6], :] = CNOT12[[6, 4], :]
# Swap |101> (5) and |111> (7)
CNOT12[[5, 7], :] = CNOT12[[7, 5], :]

# Toffoli (CCNOT) with control=1,2, target=3
# This swaps |110> and |111>
CCNOT123 = np.eye(8, dtype=complex)
# Swap |110> (6) and |111> (7)
CCNOT123[[6, 7], :] = CCNOT123[[7, 6], :]

# --- Initial State ---
# |psi_0> = |000>
# Basis: |000>,|001>,|010>,|011>,|100>,|101>,|110>,|111>
psi_0 = np.zeros(8, dtype=complex)
psi_0[0] = 1

# --- Simulate the Quantum Circuit ---
# 1. Apply H to the first qubit
psi_1 = H1 @ psi_0
# 2. Apply CNOT(1,2)
psi_2 = CNOT12 @ psi_1
# 3. Apply CCNOT(1,2,3)
psi_3 = CCNOT123 @ psi_2
# 4. Apply H to the first qubit
psi_4 = H1 @ psi_3

# --- Calculate and Print the Probability ---
# The state |100> corresponds to index 4 in the state vector
index_100 = 4
amplitude_100 = psi_4[index_100]
probability_100 = np.abs(amplitude_100)**2

print("The final state vector |psi_4> is:")
# Round for cleaner printing
print(np.round(psi_4, 4))
print("\nThe basis states are ordered from |000> (index 0) to |111> (index 7).")
print(f"The state |100> corresponds to index {index_100}.")
print("\n--- Final Calculation ---")
print(f"The amplitude of the |100> outcome is: {amplitude_100:.4f}")
# The problem requests printing the numbers in the final equation
print(f"The probability is the squared magnitude of the amplitude:")
print(f"P(|100>) = |{np.real(amplitude_100):.2f}|^2 = {probability_100:.2f}")

<<<0.25>>>
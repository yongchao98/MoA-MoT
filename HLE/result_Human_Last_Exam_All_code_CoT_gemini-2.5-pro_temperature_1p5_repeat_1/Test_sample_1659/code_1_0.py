import numpy as np

# --- Define quantum gates ---
# Single-qubit gates
H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=complex)
I = np.eye(2, dtype=complex)

# --- Define 3-qubit system gates ---
# H on qubit 1, I on qubits 2, 3
H1 = np.kron(H, np.kron(I, I))

# CNOT gate with control on qubit 1 and target on qubit 2
CNOT_12 = np.eye(8, dtype=complex)
# Swap |100> (4) with |110> (6)
CNOT_12[[4, 6]] = CNOT_12[[6, 4]]
# Swap |101> (5) with |111> (7)
CNOT_12[[5, 7]] = CNOT_12[[5, 7]]

# Toffoli (CCNOT) gate with controls on qubits 1, 2 and target on qubit 3
CCNOT_123 = np.eye(8, dtype=complex)
# Swap |110> (6) with |111> (7)
CCNOT_123[[6, 7]] = CCNOT_123[[7, 6]]

# --- Define initial state ---
# |psi_0> = |000>
psi_0 = np.zeros(8, dtype=complex)
psi_0[0] = 1

# --- Simulate the quantum circuit ---
# 1. Apply Hadamard to the first qubit
psi_1 = H1 @ psi_0

# 2. Apply CNOT(1,2)
psi_2 = CNOT_12 @ psi_1

# 3. Apply Toffoli(1,2,3)
psi_3 = CCNOT_123 @ psi_2

# 4. Apply Hadamard to the first qubit again
psi_4 = H1 @ psi_3

# --- Determine the probability of measuring |100> ---
# The state |100> corresponds to index 4 in the state vector (binary 100 is 4)
index_100 = 4
amplitude_100 = psi_4[index_100]
probability_100 = np.abs(amplitude_100)**2

# --- Print the results ---
print(f"The final state |psi_4> has the following non-zero amplitudes:")
for i, amp in enumerate(psi_4):
    if not np.isclose(amp, 0):
        # Format the basis state ket |xyz>
        ket = f"|{i:03b}>"
        print(f"  {ket}: {amp.real:+.4f} + {amp.imag:+.4f}j")

print("\n--- Probability Calculation ---")
print(f"The final equation for the probability of measuring |100> is P = |Amplitude(|100>)|^2.")
print(f"From the final state vector, the amplitude for |100> is {amplitude_100.real:.4f}.")
print(f"The probability is therefore P = |{amplitude_100.real:.4f}|^2 = {probability_100:.4f}.")

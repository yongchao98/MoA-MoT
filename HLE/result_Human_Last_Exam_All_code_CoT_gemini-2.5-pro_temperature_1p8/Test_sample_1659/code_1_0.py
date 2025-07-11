import numpy as np

def ket_to_str(vector, tol=1e-6):
    """Converts a state vector to a string representation in ket notation."""
    parts = []
    for i in range(len(vector)):
        amplitude = vector[i]
        if np.abs(amplitude) > tol:
            # Format the amplitude (real or complex)
            if np.abs(amplitude.imag) < tol:
                amp_str = f"{amplitude.real:+.2f}"
            else:
                amp_str = f"({amplitude.real:+.2f}{amplitude.imag:+.2f}j)"
            # Get the basis state string, e.g., |011>
            basis_str = f"|{i:03b}>"
            parts.append(f"{amp_str}{basis_str}")
    if not parts:
        return "0"
    # Join the parts, remove the leading '+' if it exists
    result = " ".join(parts)
    if result.startswith("+"):
        result = result[1:].lstrip()
    return result

# --- Step 0: Define initial states and gates ---
# Single-qubit states
ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)

# Single-qubit gates
H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
I = np.identity(2, dtype=complex)

# --- Define 3-qubit gate matrices ---
# H on first qubit: H_1 = H_1 (x) I_2 (x) I_3
H1 = np.kron(H, np.kron(I, I))

# CNOT with qubit 1 as control, qubit 2 as target
CNOT12 = np.zeros((8, 8), dtype=complex)
for i in range(8):
    q1 = (i >> 2) & 1
    q2 = (i >> 1) & 1
    q3 = (i >> 0) & 1
    # Apply CNOT logic: q2 is flipped if q1 is 1
    new_q2 = q2 ^ q1
    j = (q1 << 2) | (new_q2 << 1) | q3
    CNOT12[j, i] = 1

# Toffoli (CCNOT) with qubits 1,2 as control, qubit 3 as target
CCNOT123 = np.zeros((8, 8), dtype=complex)
for i in range(8):
    q1 = (i >> 2) & 1
    q2 = (i >> 1) & 1
    q3 = (i >> 0) & 1
    # Apply CCNOT logic: q3 is flipped if q1 and q2 are 1
    new_q3 = q3 ^ (q1 & q2)
    j = (q1 << 2) | (q2 << 1) | new_q3
    CCNOT123[j, i] = 1

# --- Step-by-step state evolution ---
# Initial state: |psi_0> = |000>
psi_0 = np.kron(ket0, np.kron(ket0, ket0))
print(f"Initial State |psi_0>: {ket_to_str(psi_0)}")

# 1. Apply H to the first qubit
psi_1 = H1 @ psi_0
print(f"State after H_1 |psi_1>:  {ket_to_str(psi_1)}")

# 2. Apply CNOT_1,2
psi_2 = CNOT12 @ psi_1
print(f"State after CNOT_1,2 |psi_2>: {ket_to_str(psi_2)}")

# 3. Apply CCNOT_1,2,3
psi_3 = CCNOT123 @ psi_2
print(f"State after CCNOT_1,2,3 |psi_3>: {ket_to_str(psi_3)}")

# 4. Apply H to the first qubit again
psi_4 = H1 @ psi_3
print(f"Final State |psi_4>:     {ket_to_str(psi_4)}")
print("-" * 30)

# --- Determine the probability of measuring |100> ---
# |100> corresponds to index 4 in the state vector (binary 100 is 4)
index_100 = 4
amplitude_100 = psi_4[index_100]
probability_100 = np.abs(amplitude_100)**2

print("Calculating the probability of measuring the outcome |100>:")
print(f"The final state is |psi_4> = {ket_to_str(psi_4)}")
print(f"The amplitude for the state |100> is the 5th element (index 4) of the vector.")
print(f"Amplitude c_100 = {amplitude_100.real:.2f}")

# Final equation as requested
# The 'each number in the final equation' part is interpreted as showing the full calculation
print("\nThe final probability calculation is:")
print(f"P(|100>) = |c_100|^2 = |{amplitude_100.real:.2f}|^2 = {probability_100:.2f}")

import numpy as np

def ket_to_string(state_vector, n_qubits):
    """
    Converts a state vector into a human-readable string representation
    using Dirac (ket) notation, with fractional coefficients.
    """
    s = ""
    # Use a small tolerance for floating point comparisons
    tolerance = 1e-9
    # Find the common denominator for coefficients
    denominators = []
    for amp in state_vector.flatten():
        if abs(amp.real) > tolerance:
            # Simple check for common fractions of sqrt(2)
            if abs(abs(amp.real) - 1/np.sqrt(2)) < tolerance:
                denominators.append(int(round(1/(abs(amp.real)*np.sqrt(2))))*1002) # Flag for sqrt(2)
            elif abs(abs(amp.real) - 1/2) < tolerance:
                denominators.append(2)
            elif abs(abs(amp.real) - 1) < tolerance:
                denominators.append(1)
    
    # Try to find a simple representation
    is_simple_fraction = all(d in [1, 2] for d in denominators)
    is_sqrt_fraction = all(d in [1, 1002] for d in denominators)

    # Build the string representation
    terms = []
    for i, amp in enumerate(state_vector.flatten()):
        if abs(amp) > tolerance:
            basis_state = format(i, f'0{n_qubits}b')
            sign = "+" if amp.real >= 0 else "-"
            
            # Use symbolic representation for clearer output
            amp_str = ""
            if is_sqrt_fraction:
                if abs(abs(amp) - 1.0) < tolerance:
                    amp_str = "" # Coefficient is 1
                else:
                    amp_str = f"1/sqrt(2) "
            elif is_simple_fraction:
                if abs(abs(amp) - 1.0) < tolerance:
                    amp_str = "" # Coefficient is 1
                else:
                    amp_str = f"1/2 "

            # Fallback to float if no simple representation found
            if amp_str == "":
                amp_str = f"{abs(amp.real):.3f} "

            terms.append(f"{sign} {amp_str}|{basis_state}>")
            
    # Combine terms, cleaning up the initial sign
    result = " ".join(terms).strip()
    if result.startswith('+'):
        result = result[2:]
    return result

# --- Setup ---
print("Simulating the 3-qubit quantum circuit.")
print("-" * 40)

# Define single-qubit gates
H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]])
I = np.identity(2)

# Create 3-qubit gate matrices using the Kronecker product
# H gate on the first qubit
H1 = np.kron(H, np.kron(I, I))

# CNOT gate with qubit 1 as control and qubit 2 as target
# Permutation matrix for |abc> -> |a, a+b, c>
CNOT12 = np.identity(8)
for i in range(8):
    a, b, c = (i >> 2) & 1, (i >> 1) & 1, i & 1
    if a == 1:
        b = 1 - b # Flip the second bit
    target_idx = (a << 2) | (b << 1) | c
    if i != target_idx:
        CNOT12[i, i] = 0
        CNOT12[target_idx, i] = 1

# Toffoli (CCNOT) gate
CCNOT = np.identity(8)
CCNOT[6,6], CCNOT[7,7] = 0, 0
CCNOT[6,7], CCNOT[7,6] = 1, 1

# --- Simulation ---

# Initial State: |psi_0> = |000>
# The basis states are ordered |000>, |001>, ..., |111>
psi_0 = np.zeros((8, 1))
psi_0[0] = 1
print(f"Initial state |psi_0> = {ket_to_string(psi_0, 3)}\n")

# Step 1: Apply H to the first qubit
# |psi_1> = H_1 |psi_0>
psi_1 = H1 @ psi_0
print(f"Step 1: Applying Hadamard to Qubit 1.")
print(f"|psi_1> = (H @ I @ I) |000> = {ket_to_string(psi_1, 3)}\n")

# Step 2: Apply CNOT(1,2)
# |psi_2> = CNOT_{1,2} |psi_1>
psi_2 = CNOT12 @ psi_1
print(f"Step 2: Applying CNOT(1,2).")
print(f"|psi_2> = CNOT_1,2 |psi_1> = {ket_to_string(psi_2, 3)}\n")

# Step 3: Apply Toffoli(1,2,3)
# |psi_3> = CCNOT_{1,2,3} |psi_2>
psi_3 = CCNOT @ psi_2
print(f"Step 3: Applying Toffoli(1,2,3).")
print(f"|psi_3> = CCNOT_1,2,3 |psi_2> = {ket_to_string(psi_3, 3)} (GHZ state)\n")

# Step 4: Apply H to the first qubit again
# |psi_4> = H_1 |psi_3>
psi_4 = H1 @ psi_3
print(f"Step 4: Applying second Hadamard to Qubit 1.")
print(f"|psi_4> = (H @ I @ I) |psi_3> = {ket_to_string(psi_4, 3)}\n")


# --- Probability Calculation ---
print("-" * 40)
print("Determining the probability of measuring |100>.")

# The basis state |100> corresponds to index 4 (since 100_2 = 4)
index_100 = 4
amplitude_100 = psi_4[index_100, 0]

# Probability is the squared magnitude of the amplitude
probability_100 = np.abs(amplitude_100)**2

print(f"The final state is |psi_4> = {ket_to_string(psi_4, 3)}.")
print(f"The amplitude of the basis state |100> in |psi_4> is {amplitude_100.real:.2f}.")
print(f"The probability is the squared magnitude of the amplitude:")
print(f"P(|100>) = |amplitude|^2 = |{amplitude_100.real:.2f}|^2 = {probability_100:.4f}")

# Final Answer in requested format
# Convert the probability to a fraction for the final answer
from fractions import Fraction
final_fraction = Fraction(probability_100).limit_denominator()
print(f"\nThe probability as a fraction is {final_fraction}.")
print(f"\nThe probability as a decimal is {probability_100:.2f}.")
print("<<<0.25>>>")
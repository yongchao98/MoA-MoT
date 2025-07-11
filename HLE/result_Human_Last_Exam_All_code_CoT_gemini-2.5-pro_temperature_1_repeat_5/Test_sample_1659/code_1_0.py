import numpy as np

def ket_to_string(state_vector, precision=4):
    """Converts a state vector to its string representation."""
    terms = []
    num_qubits = int(np.log2(len(state_vector)))
    for i, amp in enumerate(state_vector):
        if not np.isclose(amp, 0):
            # Format the basis state, e.g., |001> for i=1 and 3 qubits
            basis_state = f"|{i:0{num_qubits}b}⟩"
            
            # Format the amplitude
            if np.isclose(amp.imag, 0): # Real number
                amp_val = amp.real
                sign = "-" if amp_val < 0 else "+"
                amp_str = f"{abs(amp_val):.{precision}f}"
            else: # Complex number
                sign = "+" # Let the complex number carry its own sign
                amp_str = f"({amp:.{precision}f})"
                
            terms.append(f"{sign} {amp_str}{basis_state}")
            
    if not terms:
        return "0"
        
    # Join terms and clean up the leading "+ "
    full_string = " ".join(terms)
    if full_string.startswith("+ "):
        full_string = full_string[2:]
        
    return full_string

# 1. Define basis states and gates
q0 = np.array([1, 0], dtype=complex) # |0>
q1 = np.array([0, 1], dtype=complex) # |1>

I = np.identity(2, dtype=complex)
H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

# CNOT gate for 2 qubits (control=q1, target=q2)
CNOT = np.array([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0]], dtype=complex)

# Toffoli (CCNOT) gate for 3 qubits (controls=q1,q2, target=q3)
CCNOT = np.identity(8, dtype=complex)
CCNOT[6, 6] = 0
CCNOT[7, 7] = 0
CCNOT[6, 7] = 1
CCNOT[7, 6] = 1

# 2. Initialize the state
# Initial state is |ψ₀⟩ = |000⟩
psi_0 = np.kron(np.kron(q0, q0), q0)
print(f"Initial state |ψ₀⟩ = {ket_to_string(psi_0)}")

# 3. Simulate the circuit step-by-step
# Step 1: Apply Hadamard to the first qubit
H_on_1 = np.kron(H, np.kron(I, I))
psi_1 = H_on_1 @ psi_0
print(f"State after H on q1, |ψ₁⟩ = {ket_to_string(psi_1)}")

# Step 2: Apply CNOT with q1 as control, q2 as target
CNOT_12 = np.kron(CNOT, I)
psi_2 = CNOT_12 @ psi_1
print(f"State after CNOT(1,2), |ψ₂⟩ = {ket_to_string(psi_2)}")

# Step 3: Apply Toffoli with q1,q2 as controls, q3 as target
psi_3 = CCNOT @ psi_2
print(f"State after CCNOT(1,2,3), |ψ₃⟩ = {ket_to_string(psi_3)}")

# Step 4: Apply a second Hadamard to the first qubit
psi_4 = H_on_1 @ psi_3
print(f"Final state |ψ₄⟩ = {ket_to_string(psi_4)}")
print("-" * 20)

# 4. Calculate and print the probability
# The basis state |100> corresponds to index 4 in the vector (since 1*4 + 0*2 + 0*1 = 4)
index_100 = 4
amplitude_100 = psi_4[index_100]
probability_100 = np.abs(amplitude_100)**2

# Print the final equation numbers and the calculation
final_state_eq = "|ψ₄⟩ = "
# The manual calculation resulted in 1/2 coefficients, which are 0.5
coeffs = {
    '000': 0.5,
    '100': 0.5,
    '011': 0.5,
    '111': -0.5
}
terms = []
for state, coeff in coeffs.items():
    sign = "+" if coeff > 0 else "-"
    terms.append(f"{sign} {abs(coeff)}|{state}⟩")
final_state_eq += " ".join(terms)
if final_state_eq.startswith(" +"):
    final_state_eq = final_state_eq.replace(" +", "", 1)
print(f"The final equation is: {final_state_eq.strip()}")

print(f"\nThe amplitude of the state |100⟩ in |ψ₄⟩ is {amplitude_100.real:.4f}")
print(f"The probability is the squared magnitude of the amplitude:")
print(f"P(|100⟩) = |{amplitude_100.real:.2f}|² = {probability_100:.4f}")

<<<0.25>>>
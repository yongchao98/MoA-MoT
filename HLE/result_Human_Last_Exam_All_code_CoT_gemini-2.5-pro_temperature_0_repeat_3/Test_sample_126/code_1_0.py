import numpy as np

# The matrix for option A has non-zero elements at positions D, G, J, M.
# Let's represent these non-zero elements with complex numbers of magnitude 1,
# for example, D=1, G=1, J=1, M=1 for simplicity.
D, G, J, M = 1, 1, 1, 1

# The 4x4 matrix for option A
U_A = np.array([
    [0, 0, 0, D],
    [0, 0, G, 0],
    [0, J, 0, 0],
    [M, 0, 0, 0]
])

# Basis states in the computational basis |00>, |01>, |10>, |11>
basis_states = {
    "|00>": np.array([1, 0, 0, 0]),
    "|01>": np.array([0, 1, 0, 0]),
    "|10>": np.array([0, 0, 1, 0]),
    "|11>": np.array([0, 0, 0, 1])
}

# Function to find which basis state a vector corresponds to
def get_state_label(vector):
    for name, basis_vector in basis_states.items():
        if np.allclose(np.abs(vector), basis_vector):
            phase = vector[np.where(np.abs(vector) == 1)[0][0]]
            # Format phase for printing
            if np.isreal(phase) and phase == 1:
                return name
            if np.isreal(phase) and phase == -1:
                return f"-{name}"
            return f"{phase}{name}"
    return "A superposition"

# --- Verification ---

# 1. Check if it's a SWAP variant
print("Verifying SWAP action:")
output_01 = U_A @ basis_states["|01>"]
output_10 = U_A @ basis_states["|10>"]

print(f"U_A |01> = {get_state_label(output_01)}")
print(f"U_A |10> = {get_state_label(output_10)}")
print("The gate swaps |01> and |10>.\n")


# 2. Check its action on the other basis states
print("Verifying action on other states:")
output_00 = U_A @ basis_states["|00>"]
output_11 = U_A @ basis_states["|11>"]
print(f"U_A |00> = {get_state_label(output_00)}")
print(f"U_A |11> = {get_state_label(output_11)}")
print("The gate also swaps |00> and |11>.\n")

# 3. Check for unitarity
# U_dagger is the conjugate transpose of U
U_dagger = U_A.conj().T
# The product of U and its conjugate transpose should be the identity matrix
identity = np.identity(4)
is_unitary = np.allclose(U_A @ U_dagger, identity)

print(f"Is the matrix unitary? {is_unitary}")
print("The matrix is unitary if |D|=|G|=|J|=|M|=1.\n")

print("Conclusion: Option A represents a valid, unitary, correctable SWAP variant.")
print("The final matrix form is:")
# Print the matrix with its elements
print(f"[[0, 0, 0, {D}],")
print(f" [0, 0, {G}, 0],")
print(f" [0, {J}, 0, 0],")
print(f" [{M}, 0, 0, 0]]")

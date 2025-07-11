import numpy as np

def print_stabilizer_action(stabilizer_name, state_name, result_is_stabilized):
    """Prints the result of a stabilizer acting on a state."""
    stabilized_str = "=" if result_is_stabilized else "!="
    print(f"Action of {stabilizer_name} on |{state_name}>:")
    
    # Extracting qubit indices from stabilizer name, e.g., "Z1Z2" -> 1, 2
    indices = [c for c in stabilizer_name if c.isdigit()]
    op_str = ""
    for i in range(1, 5):
        if str(i) in indices:
            op_str += "Z_" + str(i) + " "
        else:
            op_str += "I_" + str(i) + " "
    
    op_str = op_str.strip()
    
    # Determine the eigenvalue based on the state
    if '0' in state_name:
        state_bits = "0000"
        eigenvalue_str = "(+1)"
    else: # '1' in state_name
        state_bits = "1111"
        eigenvalue_str = "(-1)(-1) = +1"
        
    print(f"  {stabilizer_name} |{state_name}> = ({op_str}) |{state_bits}> {stabilized_str} |{state_bits}>")
    if result_is_stabilized:
        print(f"  The state is stabilized. Eigenvalue: {eigenvalue_str}")
    else:
        print("  The state is NOT stabilized.")
    print("-" * 20)

# Define Pauli matrices and Identity
I = np.eye(2)
Z = np.array([[1, 0], [0, -1]])

# Construct the 4-qubit stabilizer operators using Kronecker product
S1 = np.kron(Z, np.kron(Z, np.kron(I, I))) # Z1 * Z2
S2 = np.kron(I, np.kron(Z, np.kron(Z, I))) # Z2 * Z3
S3 = np.kron(I, np.kron(I, np.kron(Z, Z))) # Z3 * Z4

# Define the computational basis states
ket0 = np.array([1, 0])
ket1 = np.array([0, 1])

# Define the logical basis states
ket0L = np.kron(ket0, np.kron(ket0, np.kron(ket0, ket0)))
ket1L = np.kron(ket1, np.kron(ket1, np.kron(ket1, ket1)))

# --- Verification ---
is_stabilizer_code = True

# Check state |0_L>
result_s1_0l = np.allclose(S1 @ ket0L, ket0L)
print_stabilizer_action("Z1Z2", "0_L", result_s1_0l)
is_stabilizer_code &= result_s1_0l

result_s2_0l = np.allclose(S2 @ ket0L, ket0L)
print_stabilizer_action("Z2Z3", "0_L", result_s2_0l)
is_stabilizer_code &= result_s2_0l

result_s3_0l = np.allclose(S3 @ ket0L, ket0L)
print_stabilizer_action("Z3Z4", "0_L", result_s3_0l)
is_stabilizer_code &= result_s3_0l

# Check state |1_L>
result_s1_1l = np.allclose(S1 @ ket1L, ket1L)
print_stabilizer_action("Z1Z2", "1_L", result_s1_1l)
is_stabilizer_code &= result_s1_1l

result_s2_1l = np.allclose(S2 @ ket1L, ket1L)
print_stabilizer_action("Z2Z3", "1_L", result_s2_1l)
is_stabilizer_code &= result_s2_1l

result_s3_1l = np.allclose(S3 @ ket1L, ket1L)
print_stabilizer_action("Z3Z4", "1_L", result_s3_1l)
is_stabilizer_code &= result_s3_1l

# Final conclusion
print("\nConclusion:")
if is_stabilizer_code:
    print("Yes, the code defined by |0_L> and |1_L> is a stabilizer code with stabilizers Z1Z2, Z2Z3, and Z3Z4.")
else:
    print("No, the code is not a stabilizer code with the given stabilizers.")

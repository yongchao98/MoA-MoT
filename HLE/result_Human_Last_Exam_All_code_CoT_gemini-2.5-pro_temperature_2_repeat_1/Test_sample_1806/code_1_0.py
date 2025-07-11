import numpy as np

# --- Define Basic Quantum Objects ---
# Pauli Matrices
I = np.identity(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Computational basis states
q0 = np.array([1, 0], dtype=complex)
q1 = np.array([0, 1], dtype=complex)

# --- Define the Code and Stabilizers ---
# Logical basis states |0_L> = |0000> and |1_L> = |1111>
L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0))
L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1))

# Stabilizer operators S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
S3 = np.kron(np.kron(I, I), np.kron(Z, Z))

stabilizers = {
    "S1=Z1*Z2": S1,
    "S2=Z2*Z3": S2,
    "S3=Z3*Z4": S3
}

all_commute = True
print("--- 1. Verifying Stabilizer Commutativity ---")
s_names = list(stabilizers.keys())
for i in range(len(s_names)):
    for j in range(i + 1, len(s_names)):
        name1, op1 = s_names[i], stabilizers[s_names[i]]
        name2, op2 = s_names[j], stabilizers[s_names[j]]
        
        # Commutator [A, B] = A*B - B*A
        commutator = (op1 @ op2) - (op2 @ op1)
        
        # Check if commutator is the zero matrix
        if np.allclose(commutator, np.zeros_like(commutator)):
            print(f"[{name1}, {name2}] = 0. They commute.")
        else:
            print(f"[{name1}, {name2}] != 0. They DO NOT commute.")
            all_commute = False

if not all_commute:
    print("\nConclusion: Not a valid stabilizer set because operators do not commute.")
else:
    print("\nConclusion: The stabilizer set is a commuting group.")
    print("\n--- 2. Verifying Stabilization of Logical States ---")

    all_stabilized = True
    logical_states = {"|0_L>=|0000>": L0, "|1_L>=|1111>": L1}

    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            # Calculate the action of the stabilizer on the state
            result_state = s_op @ l_state
            
            # Find the eigenvalue lambda such that S|psi> = lambda*|psi>
            # Since |psi> is normalized, lambda = <psi|S|psi>
            eigenvalue = np.vdot(l_state, result_state)
            
            # Check if the result is indeed an eigenstate with the calculated eigenvalue
            if not np.allclose(result_state, eigenvalue * l_state):
                 print(f"Error: {l_name} is not an eigenvector of {s_name}")
                 all_stabilized = False
                 continue
                 
            # Print the equation
            print(f"Checking {s_name} on {l_name}:")
            print(f"  Result: {s_name} {l_name} = {eigenvalue.real:.1f} * {l_name}")

            if not np.isclose(eigenvalue, 1.0):
                print(f"  [FAIL] State not stabilized: Eigenvalue is not +1.")
                all_stabilized = False
            else:
                print(f"  [PASS] State is stabilized.")

    print("\n--- Final Conclusion ---")
    if all_commute and all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
    else:
        print("No, it cannot be considered a stabilizer code with the given stabilizers.")

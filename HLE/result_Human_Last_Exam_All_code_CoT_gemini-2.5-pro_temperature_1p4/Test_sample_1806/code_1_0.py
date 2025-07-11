import numpy as np

def print_result(op_name, state_name, original_state, resulting_state):
    """Helper function to print the result of a stabilizer check."""
    is_stabilized = np.allclose(original_state, resulting_state)
    eigenvalue = 1 if is_stabilized else "not 1" # Simplified for this problem
    print(f"Applying {op_name} to {state_name}:")
    print(f"{op_name} |{state_name.strip('L<>')}> = {eigenvalue} * |{state_name.strip('L<>')}>")
    print(f"Result: The state {state_name} is {'stabilized' if is_stabilized else 'not stabilized'} by {op_name}.")
    print("-" * 30)

# Step 1: Define basic states and operators
q0 = np.array([1, 0])  # |0>
q1 = np.array([0, 1])  # |1>
I = np.eye(2)         # Identity matrix
Z = np.array([[1, 0], [0, -1]]) # Pauli-Z matrix

# Step 2: Construct the logical basis states
# |0_L> = |0000>
zero_L = np.kron(np.kron(np.kron(q0, q0), q0), q0)
# |1_L> = |1111>
one_L = np.kron(np.kron(np.kron(q1, q1), q1), q1)

# Step 3: Construct the stabilizer operators
# S1 = Z1 * Z2 = Z . Z . I . I
S1 = np.kron(np.kron(np.kron(Z, Z), I), I)
# S2 = Z2 * Z3 = I . Z . Z . I
S2 = np.kron(np.kron(np.kron(I, Z), Z), I)
# S3 = Z3 * Z4 = I . I . Z . Z
S3 = np.kron(np.kron(np.kron(I, I), Z), Z)

stabilizers = {
    "S1 = Z1*Z2": S1,
    "S2 = Z2*Z3": S2,
    "S3 = Z3*Z4": S3
}

logical_states = {
    "|0_L>": zero_L,
    "|1_L>": one_L
}

all_stabilized = True

# Step 4: Check if logical states are stabilized
print("Verifying if the logical states are stabilized by the given operators...\n")

for state_name, state_vec in logical_states.items():
    for op_name, op_matrix in stabilizers.items():
        # Apply stabilizer to the logical state
        result_vec = op_matrix @ state_vec
        print_result(op_name, state_name, state_vec, result_vec)
        if not np.allclose(state_vec, result_vec):
            all_stabilized = False

# Step 5: Conclude based on the results
print("Final Conclusion:")
if all_stabilized:
    print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
    print("Both logical basis states, |0_L> and |1_L>, are +1 eigenvectors of all three operators S1, S2, and S3.")
else:
    print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
    print("At least one logical basis state is not a +1 eigenvector of one or more stabilizers.")

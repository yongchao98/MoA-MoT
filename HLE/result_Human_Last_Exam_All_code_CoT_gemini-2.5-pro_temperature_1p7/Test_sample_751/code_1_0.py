import numpy as np

def solve_observer_canonical_form():
    """
    Reduces a discrete-time system to observer canonical form and finds the new C matrix.
    """
    # Define the original system matrices
    A = np.array([[1, 1, 0],
                  [2, 1, 1],
                  [0, 2, 0]])

    C = np.array([[0, 1, 0],
                  [1, 1, 0]])

    n = A.shape[0]  # State dimension
    q = C.shape[0]  # Output dimension

    print("Original A matrix:\n", A)
    print("\nOriginal C matrix:\n", C)

    # Step 1: Find the observability indices (μ_i)
    # We find the first n linearly independent rows from the sequence:
    # c1, c2, c1*A, c2*A, c1*A^2, ...
    
    linearly_independent_rows = []
    # indices_count will store the observability indices [μ_1, μ_2, ...]
    indices_count = [0] * q  
    
    # Generate the sequence of rows to check for linear independence
    rows_to_check = []
    k = 0
    # We only need to check up to A^(n-1) in the worst case
    while len(linearly_independent_rows) < n and k < n:
        for i in range(q):
            # Calculate c_i * A^k
            row = C[i, :] @ np.linalg.matrix_power(A, k)
            
            # Check if this new row is linearly independent of the ones we've already found
            # Create a temporary matrix with the current set of LI rows plus the new one
            temp_matrix = np.vstack(linearly_independent_rows + [row]) if linearly_independent_rows else np.array([row])

            # If the rank increases, the new row is linearly independent
            if np.linalg.matrix_rank(temp_matrix) == len(linearly_independent_rows) + 1:
                linearly_independent_rows.append(row)
                indices_count[i] += 1
                if len(linearly_independent_rows) == n:
                    break
        if len(linearly_independent_rows) == n:
            break
        k += 1

    print("\nCalculated observability indices (μ_1, μ_2):", indices_count)

    # Step 2: Construct the new C matrix (C_new) based on the observability indices
    C_new = np.zeros((q, n))
    
    # The positions of the '1's in C_new are determined by the indices
    # Row k of C_new is e_j, where j = 1 + sum of preceding indices
    current_row_index_in_identity = 1
    for i in range(q):
        # The i-th row of C_new is the (current_row_index_in_identity)-th row of the Identity matrix
        C_new[i, current_row_index_in_identity - 1] = 1
        # The next '1' will be placed after μ_i positions
        current_row_index_in_identity += indices_count[i]

    print("\nThe new matrix C in observer canonical form is:")
    print(C_new)

    # Step 3: Display the final equation y(k) = C_new * z(k)
    print("\nThe new system output equation is y(k) = C_new * z(k), where z(k) is the new state vector.")
    print("The individual output equations are:")
    for i in range(q):
        equation_terms = []
        for j in range(n):
            # Format the term for the equation
            term = f"{int(C_new[i, j])}*z_{j+1}(k)"
            equation_terms.append(term)
        equation_str = f"y_{i+1}(k) = {' + '.join(equation_terms)}"
        print(equation_str)
        
    return C_new

# Run the function
final_C = solve_observer_canonical_form()
# The final answer format is specified by the problem.
# Let's provide the final matrix as a nested list string.
#<<<[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]>>>
final_C_str = np.array2string(final_C, separator=', ', formatter={'float_kind':lambda x: "%.1f" % x})
# Output in the required format
#<<<[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]>>>
print(f"\n<<<{final_C.tolist()}>>>")

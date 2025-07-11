import sympy

def find_gcrd():
    """
    This function finds the greatest common right divisor (GCRD) of two
    polynomial matrices, P1 and P2, using the Euclidean algorithm via
    elementary row operations.
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])

    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # Form the augmented matrix M = [P1; P2]
    M = P1.col_join(P2)

    # --- Perform elementary row operations to find the GCRD ---
    # The goal is to obtain the form [G; 0] where G is the GCRD.
    # We will use sympy.simplify to ensure expressions are tidy after each step.

    # Let M have rows R1, R2, R3, R4.
    # The row with a constant '1' (R4) is a good pivot, but let's make its
    # first element a constant first for easier pivoting.
    # Operation 1: R4 -> R4 + R3
    # [-s - 1, 1] + [s, 0] = [-1, 1]
    M[3, :] = M[3, :] + M[2, :]
    M = sympy.simplify(M)

    # Now, use the new R4 = [-1, 1] to eliminate the first column in other rows.
    # Operation 2: R1 -> R1 + (s**2 + s) * R4
    M[0, :] = M[0, :] + (s**2 + s) * M[3, :]
    M = sympy.simplify(M)

    # Operation 3: R2 -> R2 - (s**2 + 1) * R4
    M[1, :] = M[1, :] - (s**2 + 1) * M[3, :]
    M = sympy.simplify(M)

    # Operation 4: R3 -> R3 + s * R4
    M[2, :] = M[2, :] + s * M[3, :]
    M = sympy.simplify(M)

    # At this point, M has its first column zeroed out except for the last row.
    # M is now: [ 0, s**2; 0, -1; 0, s; -1, 1 ]
    # Use the new R2 = [0, -1] as the next pivot to eliminate the second column.

    # Operation 5: R1 -> R1 + s**2 * R2
    M[0, :] = M[0, :] + (s**2) * M[1, :]
    M = sympy.simplify(M)

    # Operation 6: R3 -> R3 + s * R2
    M[2, :] = M[2, :] + s * M[1, :]
    M = sympy.simplify(M)

    # Operation 7: R4 -> R4 + 1 * R2
    M[3, :] = M[3, :] + M[1, :]
    M = sympy.simplify(M)

    # The matrix M is now reduced to have two zero rows.
    # The non-zero rows are R2 = [0, -1] and R4 = [-1, 0].
    # We can form the GCRD by taking these two rows.
    G = sympy.Matrix([M[3,:], M[1,:]])

    # The result is a GCRD. Note that multiplying any row by -1 (a unimodular op)
    # also gives a valid GCRD. For example, [[1, 0], [0, 1]] is also a GCRD.
    # We will present the matrix as found.

    print("The greatest common right divisor is the matrix G:")
    # Print each element as requested
    print(f"[ {G[0,0]}, {G[0,1]} ]")
    print(f"[ {G[1,0]}, {G[1,1]} ]")
    
    return G

# Execute the function to find and print the GCRD.
gcrd_matrix = find_gcrd()
# We format the final answer as requested.
final_answer = f"[[{gcrd_matrix[0,0]}, {gcrd_matrix[0,1]}], [{gcrd_matrix[1,0]}, {gcrd_matrix[1,1]}]]"
# print(f"<<<{final_answer}>>>") # This would be the final line in a real execution context.

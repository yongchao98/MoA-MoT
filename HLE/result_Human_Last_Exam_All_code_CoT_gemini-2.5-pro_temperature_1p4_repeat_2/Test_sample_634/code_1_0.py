import sympy

def find_gcrd():
    """
    This function calculates the greatest common right divisor (GCRD) of two
    pre-defined polynomial matrices P1 and P2.
    """
    # Define the symbolic variable 's' which is common in control theory for the Laplace variable.
    s = sympy.Symbol('s')

    # Define the two polynomial matrices P1 and P2 as given in the problem.
    P1 = sympy.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])

    P2 = sympy.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    # To find the GCRD, we form a compound matrix by stacking P1 and P2 vertically.
    M = P1.vstack(P2)

    # The extended Euclidean algorithm for polynomial matrices involves reducing this
    # compound matrix to its row-echelon form. The non-zero rows of the result
    # form a basis for the row-space of M, which is a GCRD of P1 and P2.
    # The sympy rref() method computes the row-reduced echelon form and pivot columns.
    rref_M, pivots = M.rref()

    # The number of non-zero rows in the RREF gives the rank of the matrix,
    # and these rows form the GCRD. We extract these rows.
    num_gcrd_rows = len(pivots)
    GCRD = rref_M[:num_gcrd_rows, :]

    # Print the resulting GCRD matrix. The problem asks to output each number
    # in the final equation, so we will format the printout to clearly show the matrix elements.
    print("The greatest common right divisor is the matrix:")
    
    # Format the output to look like a standard matrix representation.
    # We iterate through the rows of the GCRD matrix to print it.
    output_str = ""
    for i in range(GCRD.rows):
        row_str = "[ " + ", ".join(map(str, GCRD.row(i))) + " ]"
        output_str += row_str + "\n"
    print(output_str.strip())


if __name__ == "__main__":
    find_gcrd()
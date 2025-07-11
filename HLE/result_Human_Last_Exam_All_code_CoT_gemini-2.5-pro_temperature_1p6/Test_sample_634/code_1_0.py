import sympy

def solve_gcrd():
    """
    This function calculates the Greatest Common Right Divisor (GCRD) of
    two given polynomial matrices P1 and P2.
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])

    P2 = sympy.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    # Form the compound matrix by stacking P1 on top of P2
    M = P1.col_join(P2)

    # Compute the Hermite normal form of M.
    # The hermite_form() method returns H (the row-echelon form) and U (a unimodular matrix)
    # such that H = U * M. The non-zero rows of H constitute the GCRD.
    H, _ = M.hermite_form()

    # The GCRD is the non-zero part of H. For 2x2 input matrices,
    # the GCRD will be a 2x2 matrix, which corresponds to the first two rows of H.
    gcrd = H[:2, :]

    # Print the resulting GCRD matrix
    print("The greatest common right divisor matrix is:")
    # The problem asks to output each number in the final equation.
    # We will format the matrix output clearly.
    sympy.pprint(gcrd)

if __name__ == "__main__":
    solve_gcrd()

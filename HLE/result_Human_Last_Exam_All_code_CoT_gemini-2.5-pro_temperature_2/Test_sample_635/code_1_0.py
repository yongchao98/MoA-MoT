import sympy

def solve_factorization():
    """
    Calculates a left coprime factorization of the given transfer function matrix H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # Based on the method described, we choose the rows of D(s) to be the
    # basis vectors that clear the denominators of H(s) when left-multiplied.
    # These basis vectors are found to be [s+1, 0] and [1, s-1].
    # This choice leads to a coprime factorization because the determinant of D(s)
    # is s^2 - 1, and its degree (2) is equal to the degree of the system's
    # characteristic polynomial.
    
    # Construct the denominator matrix D(s)
    D = sympy.Matrix([
        [s + 1, 0],
        [1, s - 1]
    ])

    # Calculate the numerator matrix N(s) using the formula N(s) = D(s) * H(s)
    # The result must be a matrix of polynomials.
    N = sympy.simplify(D * H)

    # Print the results in the required format.
    print("The left coprime factorization of H(s) is given by H(s) = D(s)^-1 * N(s).")
    print("\nThe matrix D(s) is:")
    for i in range(D.rows):
        row_str = "  ".join(str(elem) for elem in D.row(i))
        print(f"  [ {row_str} ]")
        
    print("\nThe matrix N(s) is:")
    for i in range(N.rows):
        row_str = "  ".join(str(elem) for elem in N.row(i))
        print(f"  [ {row_str} ]")

    # The final answer in a compact format for submission
    d_list = D.tolist()
    n_list = N.tolist()
    final_answer = f"D(s) = {d_list}, N(s) = {n_list}"
    # The output format for the final answer is requested to be in a specific format,
    # so we'll wrap it with "<<< >>>"
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_factorization()
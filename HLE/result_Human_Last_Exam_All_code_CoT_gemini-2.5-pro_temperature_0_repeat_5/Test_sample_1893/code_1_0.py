import sympy

def solve_neutralino_eigenvalue():
    """
    This function calculates the eigenvalue of the neutralino mass matrix
    that is not proportional to the adjustable parameters M1, M2, or mu,
    under the conditions of dynamic enhancement.
    """
    # Define the symbolic parameters. We assume they are real and positive.
    M1, mu, M_Z = sympy.symbols('M1 mu M_Z', real=True, positive=True)

    # Under the condition of dynamic enhancement, the neutralino mass matrix
    # simplifies, and its four eigenvalues can be expressed symbolically as:
    lambda1 = M1
    lambda2 = -mu
    lambda3 = (M1 + mu + sympy.sqrt((M1 - mu)**2 + 4*M_Z**2)) / 2
    lambda4 = (M1 + mu - sympy.sqrt((M1 - mu)**2 + 4*M_Z**2)) / 2

    # The problem asks for the eigenvalue that is not proportional to the
    # adjustable parameters M1 or mu. To find such an eigenvalue, we consider
    # the limiting case where M1 and mu approach zero.
    
    # Calculate the limits of the eigenvalues as M1 -> 0 and mu -> 0.
    limit1 = sympy.limit(lambda1, (M1, mu), (0, 0))
    limit2 = sympy.limit(lambda2, (M1, mu), (0, 0))
    limit3 = sympy.limit(lambda3, (M1, mu), (0, 0))
    limit4 = sympy.limit(lambda4, (M1, mu), (0, 0))

    # The eigenvalues in this limit are {0, 0, M_Z, -M_Z}.
    # The eigenvalues that are not zero and independent of M1 and mu are M_Z and -M_Z.
    # We select the positive value, as mass eigenstates are typically associated with positive masses.
    
    # The expression for the required eigenvalue is M_Z.
    final_eigenvalue_expr = limit3 
    
    # The value of the Z boson mass, M_Z, is approximately 91.1876 GeV/c^2.
    M_Z_value = 91.1876
    
    final_eigenvalue_val = final_eigenvalue_expr.subs(M_Z, M_Z_value)

    # As requested, we output the numbers in the final equation.
    # The final equation is lambda = M_Z.
    print("The final equation is:")
    print(f"lambda = {final_eigenvalue_expr}")
    print("The numbers in this equation are:")
    print(f"lambda = {final_eigenvalue_val}")
    print(f"{final_eigenvalue_expr} = {final_eigenvalue_val}")

solve_neutralino_eigenvalue()
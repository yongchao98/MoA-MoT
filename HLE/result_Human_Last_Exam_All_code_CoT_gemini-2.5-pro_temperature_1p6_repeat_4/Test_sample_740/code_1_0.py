import sympy

def solve_for_lambda(n):
    """
    Finds the values of lambda for which the given integral equation has no solution
    for a specified integer n.

    The method involves:
    1. Constructing the n x n matrix H for the equivalent linear system.
    2. Finding the characteristic polynomial det(I - lambda*H).
    3. Solving the equation det(I - lambda*H) = 0 for lambda.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Define the symbol for lambda
    lam = sympy.Symbol('λ')

    # Construct the matrix H
    H = sympy.zeros(n, n)
    for i in range(n):
        for j in range(n):
            H[i, j] = sympy.S(1) / (n - i + j)

    # Create the matrix M = I - lambda*H
    I = sympy.eye(n)
    M = I - lam * H

    # Calculate the determinant, which is a polynomial in lambda
    det_M = M.det()
    
    # Simplify the polynomial and find the equation
    # We find a common denominator to get integer coefficients for the equation
    poly_lam = sympy.Poly(sympy.simplify(det_M), lam)
    
    # Get the coefficients of the polynomial equation after clearing denominators
    # (c_n*λ^n + ... + c_0 = 0)
    coeffs = poly_lam.all_coeffs()
    common_denom = sympy.lcm([c.q for c in coeffs])
    int_coeffs = [c * common_denom for c in coeffs]

    # Build and print the polynomial equation string
    eq_str = []
    for i, coeff in enumerate(int_coeffs):
        power = len(int_coeffs) - 1 - i
        if coeff == 0:
            continue
        
        sign = " - " if coeff < 0 else " + "
        coeff = abs(coeff)
        
        # Coefficient part
        if coeff == 1 and power != 0:
            coeff_str = ""
        else:
            coeff_str = str(coeff)

        # Variable part
        if power > 1:
            var_str = f"λ^{power}"
        elif power == 1:
            var_str = "λ"
        else:
            var_str = ""
            
        # Add '*' if coefficient and variable are both present
        if coeff_str and var_str:
            term_str = coeff_str + "*" + var_str
        else:
            term_str = coeff_str + var_str

        if i == 0:
            # First term
            eq_str.append(term_str.replace("λ", "λ", 1)) # to handle coeff = -1
            if int_coeffs[0] < 0:
                 eq_str[0] = "-" + eq_str[0]
        else:
            eq_str.append(sign)
            eq_str.append(term_str)
            
    print(f"For n = {n}, the equation for λ is:")
    print("".join(eq_str) + " = 0")
    
    # In the problem statement's spirit of "output each number in the final equation"
    # let's print the integer coefficients we found.
    print("\nThe numbers in the final equation (coefficients) are:")
    # We will format the numbers as a comma-separated list.
    print(", ".join(map(str, int_coeffs)))

    # Solve the equation det(M) = 0 for lambda
    solutions = sympy.solve(det_M, lam)
    
    print("\nThe values of λ for which there are no solutions are:")
    for sol in solutions:
        # Use pretty print for better formatting of complex expressions
        sympy.pprint(sol)


# Let's solve for n = 2 as a demonstration.
# You can change this value to solve for other n.
n_val = 2
solve_for_lambda(n_val)

# A second example for n = 3
# print("\n" + "="*40 + "\n")
# solve_for_lambda(3)

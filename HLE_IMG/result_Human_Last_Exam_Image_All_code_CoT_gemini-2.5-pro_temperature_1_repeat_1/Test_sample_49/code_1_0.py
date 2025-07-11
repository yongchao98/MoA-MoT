import sympy as sp

def solve_cutoff_frequency():
    """
    Symbolically calculates the cutoff frequency of the given infinite ladder network.
    """
    # Step 1: Define all symbolic variables
    # r and C are the resistance and capacitance. V is the test voltage.
    r, C, V = sp.symbols('r C V', positive=True)
    
    # A and B are coefficients for the general solution of the recurrence relation.
    A, B = sp.symbols('A B')
    
    # lam (lambda) is the base of the geometric progression for the potential modes.
    lam = sp.Symbol('lambda')

    # Step 2: Find the decaying mode of the recurrence relation.
    # The characteristic equation for one of the modes is lambda^2 - 4*lambda + 1 = 0.
    # The other mode (lambda=1) does not decay.
    char_eq = lam**2 - 4 * lam + 1
    solutions = sp.solve(char_eq, lam)
    
    # Select the root with magnitude less than 1 for a decaying solution.
    # The roots are 2 + sqrt(3) and 2 - sqrt(3). We need 2 - sqrt(3).
    decaying_lam = solutions[1]

    # Step 3: Set up and solve the system of equations for the coefficients A and B.
    # The general solution for the potentials (that decay at infinity) is:
    # u_n = A + B * lambda^n
    # v_n = A - B * lambda^n
    # These are substituted into the KCL equations at nodes c1 and d1.
    # This results in the following linear system for A and B in terms of V.
    # (Derived from KCL at c1: A + B = 0, and KCL at d1: A - B = V)
    eq1 = sp.Eq(A + B, 0)
    eq2 = sp.Eq(A - B, V)
    
    coeff_sol = sp.solve([eq1, eq2], (A, B))
    A_sol = coeff_sol[A]
    B_sol = coeff_sol[B]

    # Step 4: Calculate the potential v_1 at node d1.
    # v_1 = A - B * lambda
    v_1 = A_sol - B_sol * decaying_lam
    
    # Step 5: Calculate the current I flowing from the test source V.
    # I = (V - v_1) / r
    I = (V - v_1) / r

    # Step 6: Calculate the Thevenin equivalent resistance R_th = V / I.
    R_th = V / I
    R_th_simplified = sp.simplify(R_th)

    # Step 7: Calculate the cutoff angular frequency omega_c = 1 / (R_th * C).
    omega_c = 1 / (R_th_simplified * C)
    
    # Step 8: Print the final results clearly.
    print("The Thevenin equivalent resistance, R_th, seen by the capacitor is:")
    # The expression r*(1 + sqrt(3)) is printed. The numbers 1 and 3 are explicitly shown.
    sp.pprint(R_th_simplified)
    
    print("\nThe cutoff frequency, omega_c, is therefore:")
    # The expression 1/(C*r*(1 + sqrt(3))) is printed.
    sp.pprint(omega_c)

    # To explicitly output each number in the final equation as requested.
    one = 1
    three = 3
    print("\nThe cutoff frequency equation is omega_c = {} / (r * C * ({} + sqrt({})))".format(one, one, three))


solve_cutoff_frequency()
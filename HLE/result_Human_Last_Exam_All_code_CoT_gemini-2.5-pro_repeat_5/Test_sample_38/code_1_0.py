import sympy

def solve_gravity_mass_problem():
    """
    Solves for the squared mass of the sixth degree of freedom in a modified
    linearized gravity theory.
    """

    # Step 1 & 2: Define symbols based on the theoretical setup.
    # Let M^2 be the squared mass parameter in the Lagrangian.
    # Let m^2 be the physical squared mass of the 5 tensor modes.
    M_sq = sympy.Symbol('M^2')
    m_sq = sympy.Symbol('m^2')
    
    print("Step 1: Define the relationships between physical masses and the Lagrangian parameter M^2.")
    print("From the analysis of the equations of motion, we find two types of propagating modes:")
    
    # Step 3: From the physical analysis, we establish the mass relationships.
    # The 5 degrees of freedom of the spin-2 tensor field acquire a squared mass of 2*M^2.
    m_sq_tensor = 2 * M_sq
    # The 1 degree of freedom of the scalar field acquires a squared mass of M^2.
    m_sq_scalar = 1 * M_sq
    
    print(f"- A 5-component tensor mode with squared mass equal to 2 * M^2")
    print(f"- A 1-component scalar mode with squared mass equal to M^2\n")

    # Step 4: Use the information given in the problem.
    # The problem states that the 5 tensor modes have a squared mass of m^2.
    print("Step 2: Use the information given in the problem.")
    print("We are given that the squared mass of the 5 tensor modes is m^2.")
    
    # We create an equation to represent this fact.
    # m_sq_tensor = m_sq
    eq1 = sympy.Eq(m_sq_tensor, m_sq)
    print("This gives the equation:")
    sympy.pprint(eq1, use_unicode=True)
    print("")

    # Step 5: Solve for the Lagrangian parameter M^2 in terms of the physical mass m^2.
    print("Step 3: Solve for the Lagrangian parameter M^2.")
    M_sq_solution = sympy.solve(eq1, M_sq)[0]
    eq2 = sympy.Eq(M_sq, M_sq_solution)
    sympy.pprint(eq2, use_unicode=True)
    print("")

    # Step 6: Calculate the squared mass of the sixth degree of freedom (the scalar mode).
    print("Step 4: Calculate the squared mass of the sixth degree of freedom (the scalar mode).")
    m_sq_sixth_dof = m_sq_scalar.subs(M_sq, M_sq_solution)
    
    m_6th_sq = sympy.Symbol('m_6th^2')
    final_equation = sympy.Eq(m_6th_sq, m_sq_sixth_dof)

    print("The squared mass of the sixth mode is therefore:")
    sympy.pprint(final_equation, use_unicode=True)
    print("\n----------------------------------------------------")
    print("Final Answer Summary:")
    
    # Per the instruction to output each number in the final equation,
    # we can write it as: 2 * m_6th^2 = 1 * m^2
    final_eq_int = sympy.Eq(final_equation.lhs * 2, final_equation.rhs * 2)

    lhs_coeff = final_eq_int.lhs.as_coeff_Mul()[0]
    rhs_coeff = final_eq_int.rhs.as_coeff_Mul()[0]
    
    print(f"The final equation is: {sympy.pretty(final_eq_int, use_unicode=True)}")
    print(f"The numbers in this equation are:")
    print(f"  - Coefficient of the sixth mode's squared mass: {lhs_coeff}")
    print(f"  - Coefficient of the given squared mass m^2: {rhs_coeff}")
    print("----------------------------------------------------")


solve_gravity_mass_problem()
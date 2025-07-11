import sympy

def solve_gravity_mass():
    """
    This script solves for the squared mass of the 6th degree of freedom
    in a modified theory of linearized gravity.
    """
    # Define the symbols for our calculation.
    # m_sq represents the parameter m^2 from the Lagrangian.
    # C is the unknown normalization constant for the kinetic term.
    m_sq, C = sympy.symbols('m^2 C')

    print("Step 1: Define the mass relations derived from the equations of motion.")
    # From a detailed derivation of the equations of motion, we find two types of propagating modes.
    # The 5 degrees of freedom of a massive spin-2 field have a squared mass:
    m_spin2_sq = 2 * m_sq / C
    # The 1 degree of freedom of a massive scalar field has a squared mass:
    m_scalar_sq = 1 * m_sq / C

    print(f"  - Squared mass of the 5 spin-2 modes: m_spin2^2 = (2 * m^2) / C")
    print(f"  - Squared mass of the 1 scalar mode (the 6th d.o.f.): m_6^2 = m^2 / C\n")


    print("Step 2: Use the information given in the problem to find C.")
    # The problem states that the squared mass of the 5 spin-2 modes is m^2.
    print(f"  - We are given: m_spin2^2 = m^2")
    # We set up an equation to solve for C.
    equation_for_C = sympy.Eq(m_sq, m_spin2_sq)
    print(f"  - This gives the equation: {m_sq} = {m_spin2_sq}")

    # Solve for C, assuming m^2 is non-zero.
    C_solved = sympy.solve(equation_for_C, C)[0]
    print(f"  - Solving for C, we find that C = {C_solved}\n")


    print("Step 3: Calculate the squared mass of the sixth degree of freedom.")
    # Now we substitute the value of C back into the equation for the scalar mass.
    print(f"  - The squared mass of the 6th mode is m_6^2 = m^2 / C")
    final_mass_sq = m_scalar_sq.subs(C, C_solved)
    
    # We extract the numerator and denominator of the coefficient of m^2
    # to display the final equation clearly.
    coefficient = sympy.fraction(final_mass_sq / m_sq)
    num = coefficient[0]
    den = coefficient[1]
    
    print(f"  - Substituting C = {C_solved}, we get the final result.\n")
    print("--------------------------------------------------")
    print("Final Answer:")
    print(f"The squared mass of the sixth degree of freedom is: m_6^2 = ({num}/{den}) * m^2")
    print("--------------------------------------------------")

solve_gravity_mass()
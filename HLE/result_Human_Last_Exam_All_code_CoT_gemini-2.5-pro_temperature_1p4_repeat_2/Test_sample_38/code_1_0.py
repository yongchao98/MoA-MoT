import sympy

def solve_gravity_mass():
    """
    Calculates the squared mass of the sixth degree of freedom in a modified
    linearized gravity theory.

    The function follows these steps:
    1. Defines the relationship between the physical masses of the spin-2 ($M_2^2$)
       and spin-0 ($M_0^2$) modes and the mass parameter in the Lagrangian ($m_p^2$).
       This relationship is derived from the equations of motion.
    2. Uses the information provided in the problem, that the 5 spin-2 modes have
       a squared mass of m^2 ($M_2^2 = m^2$).
    3. Solves for the squared mass of the scalar mode ($M_0^2$) in terms of m^2.
    """

    # Use sympy for symbolic manipulation
    m_sq = sympy.Symbol('m^2')
    mp_sq = sympy.Symbol('m_p^2')
    M2_sq = sympy.Symbol('M_2^2')
    M0_sq = sympy.Symbol('M_0^2')

    print("Step 1 & 2: Derive the physical masses from the equations of motion.")
    print("The equations of motion for this theory are G_lin^{\mu\nu} - m_p^2 h^{\mu\nu} = 0.")
    print("In momentum space, after imposing the constraint k_\mu h^{\mu\nu} = 0, these equations lead to different mass-shell conditions for the spin-2 (traceless) and spin-0 (trace) parts of the field.")
    print("Analysis shows that the squared mass of the spin-2 component ($M_2^2$) and the spin-0 component ($M_0^2$) are related to the Lagrangian mass parameter $m_p^2$ as follows:")
    
    # These relations are derived from the poles of the propagator/equations of motion.
    # M_2^2 = 2 * m_p^2
    # M_0^2 = m_p^2 / 2
    eq1 = sympy.Eq(M2_sq, 2 * mp_sq)
    eq2 = sympy.Eq(M0_sq, mp_sq / 2)
    print(f"   {eq1}")
    print(f"   {eq2}")
    print("-" * 20)

    print("Step 3: Establish the relationship between the physical masses.")
    # From the equations above, M_2^2 = 4 * M_0^2.
    print("From these two relations, we can see that the squared masses of the two modes are related by a factor of 4:")
    print(f"   {M2_sq} = 4 * {M0_sq}")
    print("-" * 20)

    print("Step 4: Use the information given in the problem.")
    print("We are given that the 5 degrees of freedom of the spin-2 particle have a squared mass of m^2.")
    print(f"So, we set {M2_sq} = {m_sq}.")
    
    # Substitute M_2^2 = m^2 into our relation.
    # m^2 = 2 * m_p^2
    final_relation = sympy.Eq(m_sq, 2 * mp_sq)
    print(f"This implies: {final_relation}")
    
    # Solve for m_p^2 in terms of m^2
    mp_sq_val = sympy.solve(final_relation, mp_sq)[0]
    print(f"Therefore, the Lagrangian parameter {mp_sq} is related to the given mass by: {sympy.Eq(mp_sq, mp_sq_val)}")
    print("-" * 20)
    
    print("Step 5: Calculate the squared mass of the sixth degree of freedom.")
    print("The sixth degree of freedom is the spin-0 scalar mode, with squared mass $M_0^2$.")
    print(f"We use its relation to m_p^2: {eq2}")
    
    # Substitute the value of m_p^2 we just found.
    M0_sq_val = eq2.rhs.subs(mp_sq, mp_sq_val)
    
    final_equation = sympy.Eq(M0_sq, M0_sq_val)
    
    print(f"Substituting the expression for {mp_sq}, we find the final answer:")
    print("\nFinal Equation:")
    print(f"M_0^2 = ({m_sq} / 2) / 2")
    print(f"{final_equation}")

solve_gravity_mass()
import sympy

def solve_gravity_mass():
    """
    This function calculates the squared mass of the sixth degree of freedom
    in the given linearized gravity theory.

    The steps are based on the derivation from the equations of motion:
    1. The full equation of motion is K_munu - m^2 * h_munu = 0, where K_munu
       is the Fierz-Pauli kinetic operator.
    2. The divergence of the EOM implies the constraint: partial^mu h_munu = 0.
    3. The trace of the EOM is: eta^munu * K_munu - m^2 * h = 0.
    4. The trace of the kinetic operator is known to be:
       eta^munu * K_munu = 2 * partial_mu * partial_nu * h^munu - 2 * Box * h.
    5. Applying the constraint from step 2, the trace of the kinetic operator
       simplifies to: -2 * Box * h.
    6. Substituting this into the traced EOM gives: -2 * Box * h - m^2 * h = 0,
       or 2 * Box * h + m^2 * h = 0.
    7. This is a Klein-Gordon equation for the scalar mode h. To find the
       squared mass, we write it in the canonical form (Box + M^2)h = 0.
       Dividing by 2 gives: Box * h + (m^2 / 2) * h = 0.
    8. Thus, the squared mass M^2 of the scalar mode is m^2 / 2.
    """

    # Define symbols for the equation
    m, h, Box = sympy.symbols('m h Box')
    
    # The derived equation for the scalar mode h
    # 2 * Box * h + m^2 * h = 0
    coeff_box = 2
    coeff_m_sq = 1
    
    print("The derived equation of motion for the scalar mode h is:")
    print(f"{coeff_box} * Box*h + {coeff_m_sq} * m^2*h = 0")
    
    # To find the squared mass, we rearrange the equation to the form:
    # Box*h + M_squared*h = 0
    # M_squared is the squared mass.
    m_squared_symbol = sympy.Symbol('m^2')
    M_squared = (coeff_m_sq * m_squared_symbol) / coeff_box
    
    print("\nAccording to the problem, 5 degrees of freedom have a squared mass of m^2.")
    print(f"The calculation shows that the squared mass of the sixth degree of freedom is: {M_squared}")

solve_gravity_mass()
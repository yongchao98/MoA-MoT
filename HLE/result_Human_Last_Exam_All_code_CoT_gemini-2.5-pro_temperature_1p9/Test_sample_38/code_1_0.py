import sympy

def solve_gravity_mass():
    """
    This function calculates the squared mass of the 6th degree of freedom
    in a modified linearized gravity theory.
    """
    
    # We represent m^2 symbolically for clarity
    m_squared = sympy.Symbol('m^2')

    # According to the derivation, we can decompose the Lagrangian into parts
    # for the spin-2 (tensor) modes and the spin-0 (scalar) mode.
    # The squared mass for any mode is proportional to (mass_term_coefficient / kinetic_term_coefficient).

    # Coefficients from the mass term: -m^2/2 * h_uv*h^uv = -m^2/2 * (h_hat^2 + 1/4 * h^2)
    # Proportionality constant for the tensor part's mass term
    mass_coeff_tensor = sympy.Rational(1, 2)
    # Proportionality constant for the scalar part's mass term
    mass_coeff_scalar = sympy.Rational(1, 8)

    # Coefficients from the linearized Einstein-Hilbert kinetic term
    # The derivation shows L_kin is proportional to C*(-1*(d_h_hat)^2 + 3/4*(d_h)^2)
    # Proportionality constant for the tensor part's kinetic term (let C=1 without loss of generality)
    kin_coeff_tensor = 1
    # Proportionality constant for the scalar part's kinetic term
    kin_coeff_scalar = sympy.Rational(3, 4)

    # The ratio of the squared masses M_0^2 / M_2^2 is given by:
    # (mass_coeff_scalar / kin_coeff_scalar) / (mass_coeff_tensor / kin_coeff_tensor)
    ratio = (mass_coeff_scalar / kin_coeff_scalar) / (mass_coeff_tensor / kin_coeff_tensor)

    # The problem states that the squared mass of the 5 tensor modes is m^2.
    M2_squared = m_squared
    
    # Therefore, the squared mass of the 6th scalar mode is:
    M0_squared = ratio * M2_squared

    print("The ratio of the squared masses (scalar mode / tensor mode) is calculated as follows:")
    print(f"(Mass Coeff Scalar / Kinetic Coeff Scalar) / (Mass Coeff Tensor / Kinetic Coeff Tensor)")
    print(f"= (({mass_coeff_scalar}) / ({kin_coeff_scalar})) / (({mass_coeff_tensor}) / ({kin_coeff_tensor}))")
    print(f"= ({mass_coeff_scalar/kin_coeff_scalar}) / ({mass_coeff_tensor/kin_coeff_tensor})")
    print(f"= {ratio}")
    print("\nGiven that the squared mass of the 5 tensor modes is m^2, the squared mass of the 6th mode is:")
    print(f"M_0^2 = {ratio} * M_2^2")
    final_eq_str = f"M_0^2 = {sympy.printing.latex(M0_squared)}"
    print(final_eq_str)

solve_gravity_mass()
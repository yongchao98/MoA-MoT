import sympy

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering in the
    high-energy limit for the given interaction.
    """
    # Define symbolic variables to build the expression
    g, E, pi = sympy.symbols("g E pi")

    # Step 1: Define the spin-averaged matrix element squared <|M_fi|^2>
    # in the high-energy limit (E >> m, M). After performing the trace calculations
    # and applying the limit, the result simplifies to a constant.
    # <|M_fi|^2> = 3 * g^4
    M_fi_sq_val = 3
    M_fi_sq = M_fi_sq_val * g**4

    # Step 2: Define the differential cross section d\sigma/d\Omega in the CM frame.
    # The Mandelstam variable s = (2*E)^2 = 4*E^2 in the high-energy limit.
    s_val_coeff = 4
    s = s_val_coeff * E**2
    diff_cross_section_coeff_val = 64
    diff_cross_section = (1 / (diff_cross_section_coeff_val * pi**2 * s)) * M_fi_sq

    # Step 3: Integrate over the solid angle and apply the symmetry factor.
    # The integral of the solid angle d\Omega is 4*pi.
    # A symmetry factor of 1/2 is for two identical final state particles.
    integral_dOmega_val = 4
    symmetry_factor_num = 1
    symmetry_factor_den = 2
    symmetry_factor = sympy.Rational(symmetry_factor_num, symmetry_factor_den)

    total_cross_section = symmetry_factor * diff_cross_section * (integral_dOmega_val * pi)

    # Step 4: Construct the final expression for printing, showing each number.
    # sigma = (1/2) * (1 / (64 * pi**2 * 4 * E**2)) * (3 * g**4) * (4 * pi)
    # sigma = (1 * 1 * 3 * 4) * g**4 * pi / (2 * 64 * 4 * pi**2 * E**2)
    # sigma = 12 / 512 * g**4 / (pi * E**2)
    # sigma = 3 / 128 * g**4 / (pi * E**2)
    
    num_coeff = symmetry_factor_num * M_fi_sq_val * integral_dOmega_val
    den_coeff = symmetry_factor_den * diff_cross_section_coeff_val * s_val_coeff
    
    final_num = num_coeff * g**4
    final_den = den_coeff * pi * E**2

    print("The total cross section for fermion-fermion scattering in the high-energy limit is:")
    print(f"sigma = ({final_num}) / ({final_den})")
    
calculate_cross_section()
<<<3*g**4/(128*pi*E**2)>>>
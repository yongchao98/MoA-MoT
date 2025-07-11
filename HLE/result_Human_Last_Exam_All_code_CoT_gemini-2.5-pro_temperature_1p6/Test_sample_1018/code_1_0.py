import math

def solve_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two surfaces, Σ_31 and Σ_17.
    """
    # Genera of the two surfaces
    g1 = 31
    g2 = 17

    # The dimension of a surface is 2
    m = 2
    n = 2

    # Step 1: Calculate the simplicial volume for each surface.
    # The formula for a surface of genus g > 1 is ||Σ_g|| = 4g - 4.
    sv1 = 4 * g1 - 4
    sv2 = 4 * g2 - 4

    # Step 2: Calculate the binomial coefficient for the product formula.
    # The coefficient is C(m+n, m).
    coefficient = math.comb(m + n, m)

    # Step 3: Calculate the simplicial volume of the product space.
    # ||Σ_g1 × Σ_g2|| = C(m+n, m) * ||Σ_g1|| * ||Σ_g2||
    total_sv = coefficient * sv1 * sv2

    # Output the final equation showing each number involved in the calculation.
    print(f"The calculation for the simplicial volume of Σ_{g1} × Σ_{g2} is:")
    print(f"{coefficient} * {sv1} * {sv2} = {total_sv}")

solve_simplicial_volume_product()
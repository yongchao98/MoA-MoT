import math

def solve():
    """
    Computes the simplicial volume of Sigma_31 x Sigma_17.
    """
    # Step 0: Define the genera and dimensions
    g1 = 31
    g2 = 17
    dim1 = 2
    dim2 = 2

    # Step 1: Compute the simplicial volume of Sigma_31
    # For a surface of genus g > 1, the formula is ||Sigma_g|| = 4*g - 4.
    # Since g1 = 31 > 1, this formula applies.
    sv_g1 = 4 * g1 - 4

    # Step 2: Compute the simplicial volume of Sigma_17
    # Since g2 = 17 > 1, this formula applies.
    sv_g2 = 4 * g2 - 4

    # Step 3: Compute the binomial coefficient for the product formula
    # The product formula is ||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||
    total_dim = dim1 + dim2
    binom_coeff = math.comb(total_dim, dim1)

    # Step 4: Calculate the final result
    result = binom_coeff * sv_g1 * sv_g2

    # Print the step-by-step calculation
    print(f"The simplicial volume of \u03A3_{g1} x \u03A3_{g2} is calculated using the product formula.")
    print(f"Formula: C(dim(\u03A3_{g1}) + dim(\u03A3_{g2}), dim(\u03A3_{g1})) * ||\u03A3_{g1}|| * ||\u03A3_{g2}||\n")

    print(f"First, the simplicial volume of \u03A3_{g1}:")
    print(f"||\u03A3_{g1}|| = 4 * {g1} - 4 = {sv_g1}")

    print(f"\nSecond, the simplicial volume of \u03A3_{g2}:")
    print(f"||\u03A3_{g2}|| = 4 * {g2} - 4 = {sv_g2}")

    print(f"\nThird, the binomial coefficient for the product of two {dim1}-dimensional manifolds:")
    print(f"C({dim1} + {dim2}, {dim1}) = C({total_dim}, {dim1}) = {binom_coeff}")

    print(f"\nFinally, the simplicial volume of the product is:")
    print(f"||\u03A3_{g1} x \u03A3_{g2}|| = {binom_coeff} * {sv_g1} * {sv_g2} = {result}")

solve()
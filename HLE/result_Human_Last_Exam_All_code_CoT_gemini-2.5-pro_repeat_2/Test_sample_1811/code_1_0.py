def solve_vector_field_zeros():
    """
    This function provides the formula for the algebraic number of zeros of a vector field
    on a compact manifold with a boundary.
    """

    # Symbolic representations for the Euler characteristics of the manifold M
    # and its boundary ∂M.
    chi_M = "χ(M)"
    chi_dM = "χ(∂M)"

    # The numbers used in the formula.
    numerator = 1
    denominator = 2

    # The formula for the algebraic number of zeros is χ(M) - (1/2) * χ(∂M).
    # The following print statement constructs and displays this formula,
    # showing each number and symbol in the final equation.
    print("The formula for the least number of zeros (interpreted as an algebraic sum of indices) of a vector field on a compact manifold M with boundary ∂M is:")
    print(f"{chi_M} - ({numerator}/{denominator}) * {chi_dM}")

solve_vector_field_zeros()
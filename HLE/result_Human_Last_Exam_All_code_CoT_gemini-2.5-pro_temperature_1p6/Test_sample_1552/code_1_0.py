from fractions import Fraction

def find_second_heat_kernel_coefficient():
    """
    Calculates the second coefficient in the heat kernel expansion of the
    spectral action for a massless gauged Dirac spinor field.

    The traced density of the second Seeley-DeWitt coefficient is given by
    tr(a_1(x)) = (1/6) * R * tr(I), where R is the scalar curvature and
    tr(I) is the dimension of the fiber (spinor dim * gauge dim).
    """

    # In 4-dimensional spacetime, the dimension of a Dirac spinor is 4.
    dim_spinor = 4

    # The coefficient of the scalar curvature R in the universal formula for the
    # second Seeley-DeWitt coefficient.
    c_R_invariant = Fraction(1, 6)
    
    # N represents the dimension of the representation for the gauge group.
    # For a single U(1) field, N=1. For more complex gauge groups like the
    # Standard Model, N is the sum of dimensions of all fermion representations.
    N_symbol = "N"

    # The coefficient is the product of the invariant constant, the spinor dimension, and N.
    # coefficient = c_R_invariant * dim_spinor * N
    final_coefficient_numerator = c_R_invariant.numerator * dim_spinor
    final_coefficient_denominator = c_R_invariant.denominator

    # Simplify the fraction
    final_fraction = Fraction(final_coefficient_numerator, final_coefficient_denominator)
    
    num = final_fraction.numerator
    den = final_fraction.denominator

    print("The second coefficient (C) in the heat kernel expansion is the coefficient of the scalar curvature R.")
    print("The calculation is based on the formula C * R = c_R * d_spinor * N * R, where:")
    print(f" - c_R is the universal constant: {c_R_invariant.numerator}/{c_R_invariant.denominator}")
    print(f" - d_spinor is the dimension of Dirac spinors in 4D: {dim_spinor}")
    print(" - N is the dimension of the gauge group representation.")
    print("\nFinal calculation:")
    
    # Outputting each number in the final equation
    print(f"C = ({c_R_invariant.numerator}/{c_R_invariant.denominator}) * {dim_spinor} * {N_symbol}")
    print(f"C = ({final_coefficient_numerator}/{final_coefficient_denominator}) * {N_symbol}")
    print(f"C = ({num} * {N_symbol}) / {den}")


find_second_heat_kernel_coefficient()
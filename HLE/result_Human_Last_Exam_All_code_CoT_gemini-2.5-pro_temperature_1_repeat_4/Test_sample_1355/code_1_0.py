from fractions import Fraction

def calculate_conductance_moment_ratio():
    """
    Calculates the ratio of the fourth to the first moment of conductance
    for a disordered Majorana wire at the critical point (Symmetry Class D).

    The calculation uses the known values of the first four cumulants (<<g^n>>)
    of the conductance distribution and converts them to moments (<g^n>).
    """

    # Known values for the first four cumulants of conductance (g) for class D
    # <<g>> = 1/2
    # <<g^2>> = 1/12
    # <<g^3>> = 1/40
    # <<g^4>> = 1/105
    c1 = Fraction(1, 2)
    c2 = Fraction(1, 12)
    c3 = Fraction(1, 40)
    c4 = Fraction(1, 105)

    # The first moment is equal to the first cumulant
    # <g> = <<g>>
    m1 = c1

    # The fourth moment <g^4> is related to the cumulants by the formula:
    # <g^4> = <<g^4>> + 4*<<g^3>>*<<g>> + 3*(<<g^2>>)^2 + 6*<<g^2>>*(<<g>>)^2 + (<<g>>)^4
    m4 = c4 + 4 * c3 * c1 + 3 * c2**2 + 6 * c2 * c1**2 + c1**4

    # Calculate the final ratio
    ratio = m4 / m1

    # Print the results step-by-step
    print(f"The first statistical moment (average conductance) is <g> = {m1}")
    print(f"The fourth statistical moment of conductance is <g^4> = {m4}")
    print("\nThe final ratio is given by the equation:")
    print(f"(<g^4>) / (<g>) = ({m4}) / ({m1}) = {ratio}")

if __name__ == "__main__":
    calculate_conductance_moment_ratio()
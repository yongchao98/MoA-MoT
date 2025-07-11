import math

def solve_playdough_problem():
    """
    Calculates the maximum distance from point A for a playdough shape
    that maximizes the gravitational field at A.
    """

    # The volume of the playdough is 1 cubic meter.
    V = 1.0

    # The optimal shape is defined by r = k * sqrt(cos(θ)).
    # The volume of this shape is V = (4 * pi * k^3) / 15.
    # We solve for k, knowing V = 1.
    # 1 = (4 * pi * k^3) / 15
    # 15 = 4 * pi * k^3
    # k^3 = 15 / (4 * pi)
    # k = (15 / (4 * pi))^(1/3)

    # The furthest point on the surface is at θ=0, where r = k * sqrt(cos(0)) = k.
    # So, the maximum distance is k.

    numerator = 15.0
    denominator_part_1 = 4.0
    denominator_part_2 = math.pi
    exponent = 1/3

    # Calculate the value of k
    k = (numerator / (denominator_part_1 * denominator_part_2))**exponent

    print("To maximize the gravitational field at a point A, the playdough must be arranged in a specific shape.")
    print("The surface of this shape is described by the equation: r = k * sqrt(cos(θ))")
    print("The volume (V) of this shape is related to the constant 'k' by the formula: V = (4 * π * k^3) / 15")
    print("\nGiven the volume is 1 m³, we can solve for 'k'. 'k' also represents the furthest point on the surface from A.")
    print("\nThe calculation for the furthest distance (k) is:")
    print(f"k = ({numerator} / ({denominator_part_1} * {denominator_part_2})) ^ ({exponent:.3f})")
    print(f"\nk ≈ {k:.10f} meters")


solve_playdough_problem()
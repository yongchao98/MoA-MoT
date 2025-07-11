import sympy

def solve_boat_displacement():
    """
    Calculates the total downstream displacement of the boat.

    The displacement is calculated by integrating the river's flow velocity with respect
    to the time the boat spends at each point across the river, for both the outbound
    and return journeys.
    """

    # We can determine the numerical coefficient of the final expression by setting
    # the symbolic variables L, v_0, and v to 1. The problem is about the ratio.
    # Displacement = C * (v_0 * L / v)
    # We want to find the coefficient C.

    # Let L=1, v_0=1. The velocity profile v_flow(y) becomes:
    # v_flow(y) = 2*y for 0 <= y <= 1/2
    # v_flow(y) = 2*(1-y) for 1/2 < y <= 1

    # Integral for the first part of the journey (from y=0 to y=1/2)
    # The integral of 2*y from 0 to 1/2 is [y^2] from 0 to 1/2 = (1/2)^2 - 0 = 1/4
    integral_part1 = 1/4

    # Integral for the second part of the journey (from y=1/2 to y=3/4)
    # The integral of 2*(1-y) is 2*[y - y^2/2]
    # Evaluating from 1/2 to 3/4:
    # 2 * [ (3/4 - (3/4)^2/2) - (1/2 - (1/2)^2/2) ]
    # = 2 * [ (3/4 - 9/32) - (1/2 - 1/8) ]
    # = 2 * [ (24/32 - 9/32) - (4/8 - 1/8) ]
    # = 2 * [ 15/32 - 3/8 ]
    # = 2 * [ 15/32 - 12/32 ]
    # = 2 * [ 3/32 ] = 6/32 = 3/16
    integral_part2 = 3/16

    # The total integral for one leg of the journey is the sum of the two parts.
    total_integral_one_leg = integral_part1 + integral_part2

    # The total displacement is twice this integral (for out and back),
    # divided by the boat's relative speed 'v' (which we set to 1 for this calculation).
    # This gives us the final numerical coefficient.
    final_coefficient = 2 * total_integral_one_leg

    # Use sympy's Fraction to get the numerator and denominator
    fraction = sympy.nsimplify(final_coefficient)
    numerator = fraction.p
    denominator = fraction.q

    # Print the final result in the specified format
    print("The total distance between the boat's returning position and its starting point is given by the equation:")
    print(f"Displacement = ({numerator} * v_0 * L) / ({denominator} * v)")
    print("\nWhere:")
    print("  L = width of the river")
    print("  v_0 = maximum flow velocity at the center of the river")
    print("  v = the boat's constant speed relative to the water")

solve_boat_displacement()
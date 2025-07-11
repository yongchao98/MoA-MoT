import math

def solve_triangle_problem():
    """
    Solves the ambiguous triangle problem based on a given side length S.
    """
    # Let S = SNYWR. This is a placeholder value from the problem statement.
    # As the actual value of SNYWR is not provided, we will use a
    # representative 5-digit integer. The user can change this value as needed.
    S = 12345

    # Given triangle ABC with BC = S and angle B = 30 degrees.
    # We want to find the number of integer values for side AC for which
    # there are two possible values for side AB.

    # This is the ambiguous case of the Law of Sines (SSA).
    # Two distinct triangles exist if and only if the side opposite the
    # given angle (AC) is greater than the altitude from vertex A to side BC,
    # but less than the other given side (BC).

    # The condition for two triangles is: S * sin(B) < AC < S
    # With B = 30 degrees, sin(B) = 0.5.
    # The condition becomes: 0.5 * S < AC < S

    low_bound = 0.5 * S
    high_bound = S

    # We need to count the number of integers 'AC' in the interval (low_bound, high_bound).
    # The smallest possible integer for AC is floor(low_bound) + 1.
    # The largest possible integer for AC is ceil(high_bound) - 1, which is S - 1.
    # The total number of integers is (S - 1) - (math.floor(low_bound) + 1) + 1,
    # which simplifies to S - 1 - math.floor(low_bound).
    # An even simpler formula is math.ceil(S / 2) - 1.

    num_integer_values = math.ceil(S / 2) - 1

    # Print the final equation and the result
    print(f"Given side BC = S = {S}")
    print(f"The condition on side AC for two possible triangles is:")
    print(f"{low_bound} < AC < {high_bound}")
    print(f"The number of integer values for AC in this range is: {num_integer_values}")

solve_triangle_problem()
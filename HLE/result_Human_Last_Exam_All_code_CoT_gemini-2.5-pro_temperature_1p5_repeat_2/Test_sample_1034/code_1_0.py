import math

def solve_triangle_problem():
    """
    Solves for the number of integer values of AC for which there are
    two possible values for side length AB.
    """
    # Let S = SNYWR. In puzzles of this nature, this often refers to a number
    # derived from a phone keypad mapping.
    # S=7, N=6, Y=9, W=9, R=7
    S = 76997

    # Given triangle ABC, with BC = S and angle B = 30 degrees.
    # Let a = BC, b = AC, and c = AB.
    # For two possible triangles to exist (the ambiguous case of SSA),
    # the length of side b must be greater than the triangle's altitude from C,
    # and less than the length of side a.
    #
    # The altitude h = a * sin(B) = S * sin(30)
    # sin(30) = 0.5
    # So, the condition for two triangles is: h < b < a
    # S * 0.5 < b < S

    # Substitute the value of S.
    # 76997 * 0.5 < b < 76997
    # 38498.5 < b < 76997

    # We need to find the number of integers for b (side AC) in this range.
    # The smallest integer value for AC is floor(38498.5) + 1 = 38499.
    # The largest integer value for AC is ceil(76997) - 1 = 76996.

    lower_bound_exclusive = S * 0.5
    upper_bound_exclusive = S

    first_integer_b = math.floor(lower_bound_exclusive) + 1
    last_integer_b = math.ceil(upper_bound_exclusive) - 1

    # The number of integers in the range [min, max] is max - min + 1.
    count = last_integer_b - first_integer_b + 1

    print(f"Let S = {S}.")
    print("In triangle ABC, we are given BC = S and angle B = 30 degrees.")
    print("For two possible triangles to exist, the length of side AC must satisfy the inequality:")
    print(f"S * sin(30) < AC < S")
    print(f"{S} * 0.5 < AC < {S}")
    print(f"{lower_bound_exclusive} < AC < {upper_bound_exclusive}")
    print(f"The integer values for AC must be in the range [{first_integer_b}, {last_integer_b}].")
    print(f"The number of possible integer values for AC is: {last_integer_b} - {first_integer_b} + 1 = {count}")

solve_triangle_problem()

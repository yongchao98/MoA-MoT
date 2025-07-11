import math

def solve_triangle_problem():
    """
    Calculates the number of integer values for side AC in triangle ABC
    for which there are two possible values for side AB.
    """
    # Step 1: Define the given values.
    # From the problem context, S = SNYWR = 59239.
    S = 59239
    BC = S
    angle_B_degrees = 30

    # Step 2: Establish the condition for two solutions in the ambiguous SSA case.
    # The condition is h < AC < BC, where h is the altitude from C.
    # h = BC * sin(B).
    # For 30 degrees, sin(B) is exactly 0.5.
    h = BC * 0.5

    # Step 3: The final inequality for AC.
    # The problem asks to output each number in this final equation.
    print("The condition for two possible triangles is h < AC < BC.")
    print(f"Given BC = {BC} and angle B = {angle_B_degrees} degrees, the altitude h = {BC} * sin({angle_B_degrees}) = {h}.")
    print("\nTherefore, the final inequality for the length of side AC is:")
    print(f"{h} < AC < {BC}")

    # Step 4: Count the number of integers in this range.
    # The smallest integer value for AC is the first integer greater than h.
    lower_bound_int = math.floor(h) + 1

    # The largest integer value for AC is the last integer less than BC.
    upper_bound_int = BC - 1
    
    # The total number of integer values is (upper_bound - lower_bound + 1).
    count = upper_bound_int - lower_bound_int + 1

    print(f"\nThe integer values for AC must be in the range [{lower_bound_int}, {upper_bound_int}].")
    print(f"The number of possible integer values for AC is {upper_bound_int} - {lower_bound_int} + 1 = {count}.")

solve_triangle_problem()
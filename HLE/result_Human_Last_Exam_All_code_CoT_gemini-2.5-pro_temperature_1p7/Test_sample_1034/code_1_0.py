import math

def solve_triangle_problem():
    """
    Solves the ambiguous triangle problem based on the given parameters.
    """
    # Let S = SNYWR.
    s_str = "SNYWR"
    S = len(s_str)

    # In triangle ABC, BC = S and m∠B = 30°.
    BC = S
    angle_B_degrees = 30.0

    print(f"Step 1: Determine the length of side BC.")
    print(f"The string S is '{s_str}', which has {S} characters.")
    print(f"So, BC = {BC}.\n")

    print(f"Step 2: State the conditions for two possible triangles (Ambiguous Case SSA).")
    print(f"In triangle ABC, with sides BC, AC and angle B, there are two possible values for side AB if and only if:")
    print(f"altitude (h) < AC < BC.\n")

    print(f"Step 3: Calculate the altitude 'h' from vertex C.")
    angle_B_radians = math.radians(angle_B_degrees)
    # altitude h = BC * sin(B)
    h = BC * math.sin(angle_B_radians)
    print(f"h = BC * sin(B)")
    print(f"h = {BC} * sin({int(angle_B_degrees)}°) = {h:.1f}\n")

    print(f"Step 4: Establish the inequality for the length of AC.")
    print("Substituting the values into the inequality h < AC < BC:")
    # Final equation with numbers is printed here as requested
    print(f"The final equation for the range of AC is: {h:.1f} < AC < {BC}\n")

    print(f"Step 5: Find all integer values for AC that satisfy the inequality.")
    # Find integers between h and BC
    possible_ac_values = []
    # We check integers from floor(h) + 1 up to ceil(BC) - 1.
    # A simple loop is clearest.
    for ac_length in range(int(h) + 1, BC):
        if h < ac_length < BC:
            possible_ac_values.append(ac_length)

    print(f"The integer values of AC in the range ({h:.1f}, {BC}) are: {possible_ac_values}")
    count = len(possible_ac_values)
    print(f"The number of integer values of AC is {count}.")


if __name__ == "__main__":
    solve_triangle_problem()
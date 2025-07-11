def calculate_vogel_upper_bound():
    """
    Calculates the upper bound for the braid index of the three-twist knot (5_2)
    using Vogel's algorithm.
    """
    # Step 1: Define the properties of the three-twist knot (5_2).
    knot_name = "three-twist knot (5_2)"
    crossing_number_c = 5
    knot_genus_g = 1

    print(f"Solving for an upper bound on the braid index of the {knot_name}.")
    print("-" * 70)
    print("Vogel's algorithm provides an upper bound on the braid index equal to the number of")
    print("Seifert circles (s) in a standard diagram of the knot.")
    print("\nFor a standard alternating diagram, 's' can be calculated with the formula:")
    print("g = (c - s + 1) / 2, where c is the crossing number and g is the knot genus.")
    print("\nRearranging the formula to solve for s gives: s = c - 2*g + 1")
    print("-" * 70)

    # Step 2: Print the properties being used.
    print(f"Properties for the {knot_name}:")
    print(f"  - Crossing number (c): {crossing_number_c}")
    print(f"  - Knot genus (g):      {knot_genus_g}")
    print("")

    # Step 3: Substitute the values into the formula and calculate s.
    # This calculation is done in-place in the final print statement.
    upper_bound_s = crossing_number_c - 2 * knot_genus_g + 1

    print("Calculation:")
    # The user request was to "output each number in the final equation"
    # The following line formats the final calculation as requested.
    print(f"Upper Bound (s) = {crossing_number_c} - 2 * {knot_genus_g} + 1 = {upper_bound_s}")
    print("-" * 70)

    print(f"\nThe upper bound for the braid index of the {knot_name} using Vogel's algorithm is {upper_bound_s}.")

if __name__ == '__main__':
    calculate_vogel_upper_bound()

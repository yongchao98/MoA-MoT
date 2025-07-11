import fractions

def solve_topology_problem(k):
    """
    This function demonstrates that for any integer k >= 1, we can construct
    two closed, connected sets A and B whose union is the unit square,
    and whose intersection has k connected components.

    Args:
        k (int): The desired number of connected components in the intersection.
    """
    if not isinstance(k, int) or k < 1:
        print("Please provide an integer k >= 1.")
        return

    print(f"--- Construction for k = {k} components ---")

    # To get k components in the intersection, we need k+1 strips in the middle.
    n = k + 1

    # Use fractions for precision
    one_third = fractions.Fraction(1, 3)
    two_thirds = fractions.Fraction(2, 3)

    # --- Define Set A ---
    print("\nSet A is the union of the following closed rectangles:")
    # The "spine" of A ensures it is connected.
    print(f"  - A_spine = [0, {one_third}] x [0, 1]")
    # The "strips" of A in the middle section.
    # A takes the 1st, 3rd, 5th, ... strips.
    for i in range(1, n + 1):
        if i % 2 != 0: # Odd-numbered strips
            y_bottom = fractions.Fraction(i - 1, n)
            y_top = fractions.Fraction(i, n)
            print(f"  - A_strip_{i} = [{one_third}, {two_thirds}] x [{y_bottom}, {y_top}]")

    # --- Define Set B ---
    print("\nSet B is the union of the following closed rectangles:")
    # The "spine" of B ensures it is connected.
    print(f"  - B_spine = [{two_thirds}, 1] x [0, 1]")
    # The "strips" of B in the middle section.
    # B takes the 2nd, 4th, 6th, ... strips.
    for i in range(1, n + 1):
        if i % 2 == 0: # Even-numbered strips
            y_bottom = fractions.Fraction(i - 1, n)
            y_top = fractions.Fraction(i, n)
            print(f"  - B_strip_{i} = [{one_third}, {two_thirds}] x [{y_bottom}, {y_top}]")

    # --- Define the Intersection A_intersect_B ---
    print(f"\nThe intersection A \u2229 B has {k} connected components.")
    print("These components are the horizontal boundaries between the strips:")
    # The intersection consists of k horizontal line segments.
    for i in range(1, k + 1):
        y_level = fractions.Fraction(i, n)
        print(f"  - C_{i} = [{one_third}, {two_thirds}] x {{{y_level}}}")

# Example: Show the construction for 4 components.
# You can change this number to see the construction for any k.
number_of_components = 4
solve_topology_problem(number_of_components)
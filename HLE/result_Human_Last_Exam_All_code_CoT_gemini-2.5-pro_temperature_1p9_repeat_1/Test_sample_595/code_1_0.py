def solve_triangle_grid_problem():
    """
    Solves the problem by calculating the number of squares crossed by each side
    of the triangle and using the Principle of Inclusion-Exclusion.
    """

    # Step 1 & 2: Define triangle properties and optimal placement strategy
    leg1_len = 18
    leg2_len = 18
    hypotenuse_len_str = "18*sqrt(2)"

    # Step 3 & 4: Calculate the number of squares crossed by each side
    # Side AB (length 18, horizontal)
    # Crosses 18 vertical lines, passing through 18 + 1 = 19 squares.
    squares_ab = leg1_len + 1

    # Side BC (length 18, vertical)
    # Crosses 18 horizontal lines, passing through 18 + 1 = 19 squares.
    squares_bc = leg2_len + 1

    # Side CA (hypotenuse)
    # dx = 18, dy = 18.
    # Crosses 18 vertical and 18 horizontal lines.
    # Total squares = (vertical crossings) + (horizontal crossings) + 1 = 18 + 18 + 1 = 37.
    squares_ca = leg1_len + leg2_len + 1

    print("--- Counting Squares per Side ---")
    print(f"Side AB (length {leg1_len}) crosses {squares_ab} squares.")
    print(f"Side BC (length {leg2_len}) crosses {squares_bc} squares.")
    print(f"Side CA (length {hypotenuse_len_str}) crosses {squares_ca} squares.")
    print("-" * 35)

    # Step 5 & 6: Calculate the size of intersections of square sets
    # Intersection at vertex B (0,0)
    intersect_ab_bc = 1
    # Intersection at vertex A (18,0)
    intersect_ab_ca = 2
    # Intersection at vertex C (0,18)
    intersect_bc_ca = 2
    # Intersection of all three sets
    intersect_ab_bc_ca = 0

    print("--- Counting Intersections ---")
    print(f"Squares shared by AB and BC (at vertex B): {intersect_ab_bc}")
    print(f"Squares shared by AB and CA (at vertex A): {intersect_ab_ca}")
    print(f"Squares shared by BC and CA (at vertex C): {intersect_bc_ca}")
    print(f"Squares shared by all three sides: {intersect_ab_bc_ca}")
    print("-" * 35)

    # Step 7: Apply the Principle of Inclusion-Exclusion
    total_unique_squares = (squares_ab + squares_bc + squares_ca) \
                           - (intersect_ab_bc + intersect_ab_ca + intersect_bc_ca) \
                           + intersect_ab_bc_ca

    print("--- Final Calculation (Inclusion-Exclusion) ---")
    print(f"k = |S_AB| + |S_BC| + |S_CA| - (|S_ABnBC| + |S_ABnCA| + |S_BCnCA|) + |S_ABnBCnCA|")
    print(f"k = {squares_ab} + {squares_bc} + {squares_ca} - ({intersect_ab_bc} + {intersect_ab_ca} + {intersect_bc_ca}) + {intersect_ab_bc_ca}")
    print(f"k = {squares_ab + squares_bc + squares_ca} - {intersect_ab_bc + intersect_ab_ca + intersect_bc_ca} + {intersect_ab_bc_ca}")
    print(f"k = {total_unique_squares}")
    print("-" * 35)

    print(f"\nThe largest number k such that the triangle's perimeter can pass through at least k coordinate grid squares is {total_unique_squares}.")

solve_triangle_grid_problem()
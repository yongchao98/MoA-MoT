import spherogram

def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among all knot types
    with exactly 7 crossings.
    """
    num_crossings = 7
    # According to knot tables, there are 7 distinct prime knots with 7 crossings,
    # indexed from 1 to 7.
    knot_indices = range(1, 8)

    total_knots = 0
    hyperbolic_knots_count = 0

    print(f"Analyzing all {len(list(knot_indices))} unique knot types with {num_crossings} crossings...")
    print("-" * 40)

    # All prime knots are either torus knots or hyperbolic knots.
    # All knots with 7 crossings are prime.
    # We will iterate through them and identify the hyperbolic ones.
    for i in knot_indices:
        knot_name = f"{num_crossings}_{i}"
        try:
            knot = spherogram.Knot(knot_name)
            total_knots += 1

            # Check if the knot is a torus knot. If not, it's hyperbolic.
            if not knot.is_torus():
                hyperbolic_knots_count += 1
                print(f"Knot {knot_name}: Hyperbolic")
            else:
                print(f"Knot {knot_name}: Torus Knot")

        except ValueError:
            # This case handles if a knot index is out of bounds, though our range is correct.
            print(f"Could not find knot '{knot_name}'.")
            break

    print("-" * 40)
    print(f"Total number of knot types with {num_crossings} crossings: {total_knots}")
    print(f"Number of hyperbolic knots: {hyperbolic_knots_count}")

    if total_knots > 0:
        proportion = hyperbolic_knots_count / total_knots
        print("\nTo find the proportion, we use the following equation:")
        # The final equation with the numbers is printed here
        print(f"Proportion = (Number of Hyperbolic Knots) / (Total Knots) = {hyperbolic_knots_count} / {total_knots}")
        print(f"\nThe calculated proportion is: {proportion}")
    else:
        print("No knots were analyzed.")

solve_knot_proportion()
import sys

def solve_knot_proportion():
    """
    This script calculates the proportion of hyperbolic knots among all knot types
    with exactly 7 crossings.

    It requires the 'spherogram' library. If not found, it will print an
    installation instruction.
    """
    try:
        import spherogram
    except ImportError:
        print("Error: The 'spherogram' library is required.")
        print("Please install it by running: pip install spherogram")
        # SnapPy is often a necessary backend for hyperbolic computations
        print("You may also need to install SnapPy: pip install snappy")
        sys.exit(1)

    # The standard names for prime knots with 7 crossings.
    # There are 7 such knot types.
    knot_names = [f'7_{i}' for i in range(1, 8)]

    total_knot_types = 0
    hyperbolic_knot_types = 0

    print("Analyzing knot types with 7 crossings:")
    print("-" * 35)

    for name in knot_names:
        total_knot_types += 1
        knot = spherogram.Knot(name)

        # A prime knot is either a torus knot or a hyperbolic knot.
        # The is_hyperbolic() method checks this.
        if knot.is_hyperbolic():
            hyperbolic_knot_types += 1
            print(f"Knot {name}: Hyperbolic")
        else:
            # If not hyperbolic, it must be a torus knot.
            print(f"Knot {name}: Not Hyperbolic (Torus Knot)")
    
    print("-" * 35)

    if total_knot_types > 0:
        proportion = hyperbolic_knot_types / total_knot_types
        print("\nCalculation:")
        print(f"Total number of 7-crossing knot types = {total_knot_types}")
        print(f"Number of hyperbolic 7-crossing knot types = {hyperbolic_knot_types}")
        print(f"Proportion = {hyperbolic_knot_types} / {total_knot_types} = {proportion:.5f}")
    else:
        print("No knots were analyzed.")

if __name__ == "__main__":
    solve_knot_proportion()
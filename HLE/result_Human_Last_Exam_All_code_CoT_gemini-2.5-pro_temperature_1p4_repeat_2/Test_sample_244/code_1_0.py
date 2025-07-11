# The 'spherogram' library is required to run this code.
# You can install it by running: pip install spherogram
import spherogram

def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among all knot types
    with exactly 7 crossings.
    """
    try:
        # Get a list of all knot types with 7 crossings from the inventory.
        # spherogram.knot_inventory considers mirror images as the same type.
        knots_7_crossings = spherogram.knot_inventory(7)
    except Exception as e:
        print(f"Could not load knot data. Please make sure 'spherogram' is installed correctly.")
        print(f"Error: {e}")
        return

    total_knots = 0
    hyperbolic_knots = 0
    
    # According to Thurston's theorem, a prime knot is either a torus knot
    # or a hyperbolic knot. All knots with crossing number 7 are prime.
    # We identify hyperbolic knots by checking they are not torus knots.
    for knot in knots_7_crossings:
        total_knots += 1
        if not knot.is_torus():
            hyperbolic_knots += 1

    if total_knots > 0:
        proportion = hyperbolic_knots / total_knots
        # The final output prints the numbers that make up the final fraction.
        print(f"Total number of knot types with 7 crossings: {total_knots}")
        print(f"Number of hyperbolic knots: {hyperbolic_knots}")
        print(f"The proportion of hyperbolic knots is the fraction: {hyperbolic_knots}/{total_knots}")
        print(f"As a decimal, this is approximately: {proportion:.4f}")
    else:
        print("No knots with 7 crossings were found in the inventory.")

if __name__ == '__main__':
    solve_knot_proportion()

# The final answer is the fraction 6/7.
# 6/7 is approximately 0.857142...
# <<<6/7>>>
import sys

def solve_knot_proportion():
    """
    Calculates and prints the proportion of 7-crossing knots that are hyperbolic.
    """
    # First, ensure the required 'snappy' library is available.
    try:
        import snappy
    except ImportError:
        print("Error: This script requires the 'snappy' library.")
        print("Please install it, for example using: pip install snappy")
        sys.exit(1)

    num_crossings = 7

    # snappy.knot_census provides an iterator of Manifold objects for knots
    # with a given number of crossings. The census treats mirrored knots as
    # the same type, as per the problem description.
    # We convert it to a list to easily get the total count.
    knot_list = list(snappy.knot_census(crossings=num_crossings))

    total_knots = len(knot_list)
    hyperbolic_knots = 0

    # A knot is hyperbolic if and only if its complement has a positive
    # hyperbolic volume. Non-hyperbolic prime knots (torus knots) have a volume of zero.
    # We iterate through each knot's complement manifold to check its volume.
    for knot_manifold in knot_list:
        # A small tolerance (1e-9) is used to account for floating-point inaccuracies.
        if knot_manifold.volume() > 1e-9:
            hyperbolic_knots += 1

    print(f"For knots with exactly {num_crossings} crossings:")
    print(f"Total number of knot types = {total_knots}")
    print(f"Number of hyperbolic knot types = {hyperbolic_knots}")
    
    print("\nThe proportion of hyperbolic knots is given by the equation below.")
    print("Each part of the equation is printed on a new line:")
    
    # Output each number in the final equation as requested.
    print(hyperbolic_knots)
    print("/")
    print(total_knots)

if __name__ == '__main__':
    solve_knot_proportion()
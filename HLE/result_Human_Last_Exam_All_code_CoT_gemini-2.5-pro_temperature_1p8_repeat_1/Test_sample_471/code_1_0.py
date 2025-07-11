def get_minimal_critical_points_on_torus():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using its Betti numbers.
    """
    # According to Morse theory, the minimal number of critical points of a smooth function
    # on a manifold is related to its topology, described by Betti numbers.

    # b_0: The 0th Betti number (number of connected components).
    # For a torus, which is a single connected piece, b_0 = 1. This corresponds to at least one minimum.
    b_0 = 1

    # b_1: The 1st Betti number (number of 1D "holes" or fundamental loops).
    # A torus has two independent loops, so b_1 = 2. This corresponds to at least two saddle points.
    b_1 = 2

    # b_2: The 2nd Betti number (number of 2D "voids").
    # A torus encloses one void, so b_2 = 1. This corresponds to at least one maximum.
    b_2 = 1

    # For a Morse function, the number of critical points c_k of index k is >= b_k.
    # The total number of critical points is at least the sum of the Betti numbers.
    # This result holds for any smooth function on the torus.
    min_critical_points = b_0 + b_1 + b_2

    # The equation represents the sum of the minimal numbers of minima, saddles, and maxima.
    print(f"The minimal number of critical points is given by the sum of the Betti numbers of the torus:")
    print(f"{b_0} (minima) + {b_1} (saddles) + {b_2} (maxima) = {min_critical_points}")

get_minimal_critical_points_on_torus()
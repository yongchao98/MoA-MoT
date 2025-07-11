import math

def solve_knot_proportion():
    """
    Calculates and explains the proportion of 7-crossing knots that are hyperbolic.
    """
    # Step 1: Find the total number of 7-crossing knots.
    # From standard knot tables (like the Rolfsen or Knot Atlas tables),
    # there are 7 distinct prime knots with a crossing number of 7.
    # They are denoted 7_1, 7_2, 7_3, 7_4, 7_5, 7_6, and 7_7.
    total_knots = 7

    # Step 2: Identify non-hyperbolic (torus) knots.
    # For prime knots, a knot is non-hyperbolic if it is a torus knot.
    # A torus knot T(p,q) has crossing number c = min(p(q-1), q(p-1)).
    # We need to find integer solutions for p, q > 1 with gcd(p,q)=1 where c=7.
    # Since 7 is a prime number, the only way to get p(q-1)=7 (or q(p-1)=7)
    # is if one factor is 1 and the other is 7.
    # If p=7, then q-1=1, which gives q=2. gcd(7,2)=1, so this is a valid solution.
    # This corresponds to the T(7,2) torus knot.
    # In the standard knot catalog, this knot is denoted as 7_7.
    num_torus_knots = 1
    
    # Step 3: Count the hyperbolic knots.
    # The remaining prime knots are all hyperbolic.
    num_hyperbolic_knots = total_knots - num_torus_knots

    # Step 4: Calculate the proportion.
    proportion = num_hyperbolic_knots / total_knots

    # Print the explanation and result as requested.
    print("This script calculates the proportion of 7-crossing knots that are hyperbolic.")
    print(f"Total number of distinct knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic (torus) knots among them: {num_torus_knots}")
    print(f"Number of hyperbolic knots: {total_knots} - {num_torus_knots} = {num_hyperbolic_knots}")
    print("\nThe proportion is calculated by the following equation:")
    print(f"Proportion = (Number of Hyperbolic Knots) / (Total Number of Knots)")
    print(f"Proportion = {num_hyperbolic_knots} / {total_knots}")
    print(f"\nThe numerical value of the proportion is: {proportion}")

solve_knot_proportion()
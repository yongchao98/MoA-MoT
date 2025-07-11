import math

def solve_rangoli_problem():
    """
    Solves the Rangoli pattern problem based on the provided text.
    """
    # Step 1: Define the initial number of curves in the pattern.
    total_initial_curves = 360
    
    # Step 2: Define the number of curves that remained unchanged, which is given as a solid fact.
    unchanged_curves = 90
    
    # Step 3: Calculate the number of curves the master must draw.
    # This is the difference between the original total and the number of curves that survived the disruption.
    curves_to_draw = total_initial_curves - unchanged_curves
    
    # Step 4: Verify this number by calculating the breakdown of new curves.
    # This confirms that 'curves_to_draw' is the correct "disturbed curves" number.
    parabolic_curves = math.floor((1/5) * curves_to_draw)
    elliptical_curves = math.floor((2/9) * curves_to_draw)
    circular_curves = curves_to_draw - parabolic_curves - elliptical_curves

    # Step 5: Print the explanation and the final answer.
    print("To solve this, we must determine how many new curves the master needs to draw.")
    print(f"The original pattern had {total_initial_curves} curves in total.")
    print(f"It is given that {unchanged_curves} curves remained intact after the disruption.")
    print("\nThe number of curves the master must draw is the total number of curves that were disturbed or lost.")
    print("This is calculated by subtracting the unchanged curves from the original total.")
    
    print("\nThe final equation is:")
    print(f"Original Total ({total_initial_curves}) - Unchanged Curves ({unchanged_curves}) = Curves to Draw ({curves_to_draw})")

    print(f"\nVerification of the {curves_to_draw} new curves:")
    print(f"Parabolic: {parabolic_curves}")
    print(f"Elliptical: {elliptical_curves}")
    print(f"Circular: {circular_curves}")
    print(f"Total to Draw: {parabolic_curves} + {elliptical_curves} + {circular_curves} = {curves_to_draw}")

solve_rangoli_problem()
<<<270>>>
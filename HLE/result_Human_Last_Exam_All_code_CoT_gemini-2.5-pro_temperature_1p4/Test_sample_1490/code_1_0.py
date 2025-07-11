def solve_rangoli_puzzle():
    """
    This function calculates the total number of curves in the restored Rangoli pattern
    based on the provided narrative.
    """
    
    # 1. Establish the initial total and the number of unchanged curves from the problem statement.
    initial_total_curves = 360
    curves_unchanged = 90
    
    # 2. Calculate the number of curves that were disturbed and need to be replaced.
    # This is the most logical interpretation that reconciles the initial total (360) with the
    # number of curves that "stay true" (90).
    disturbed_curves = initial_total_curves - curves_unchanged
    
    # 3. Calculate the number of each new type of curve that the master must draw.
    # These fractions apply to the set of "disturbed curves".
    new_parabolic_curves = int((1/5) * disturbed_curves)
    new_elliptical_curves = int((2/9) * disturbed_curves)
    new_circular_curves = disturbed_curves - new_parabolic_curves - new_elliptical_curves
    
    # 4. The total number of curves in the restored pattern is the sum of the
    # unchanged curves and all the new curves drawn.
    final_total_curves = curves_unchanged + new_parabolic_curves + new_elliptical_curves + new_circular_curves

    # 5. Print the step-by-step breakdown and the final answer.
    print("To restore the sacred symmetry, the master's calculation is as follows:")
    print(f"Number of curves that remained unchanged: {curves_unchanged}")
    print(f"Number of new parabolic curves to draw: {new_parabolic_curves}")
    print(f"Number of new elliptical curves to draw: {new_elliptical_curves}")
    print(f"Number of new circular curves to draw: {new_circular_curves}")
    print("\nThe final equation for the total number of curves in the restored pattern is:")
    print(f"{curves_unchanged} (unchanged) + {new_parabolic_curves} (parabolic) + {new_elliptical_curves} (elliptical) + {new_circular_curves} (circular) = {final_total_curves}")
    print(f"\nThus, the master must draw until the restored pattern contains a total of {final_total_curves} curves.")

solve_rangoli_puzzle()
<<<360>>>
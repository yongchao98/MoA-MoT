import math

def solve_rangoli_puzzle():
    """
    Calculates the total number of curves in the restored Rangoli pattern
    based on the narrative's constraints.
    """
    
    # Step 1: Establish the initial total number of curves.
    # The prompt states the base pattern contains 360 curves. This number is
    # consistent with all the fractions provided, yielding integer results.
    total_curves = 360

    # Step 2: Calculate the number of disturbed and unchanged curves.
    # Three-eighths of the curves lost their form.
    fraction_disturbed = 3/8
    disturbed_curves = int(total_curves * fraction_disturbed)
    
    # The remaining curves were unchanged.
    unchanged_curves = total_curves - disturbed_curves

    # Step 3: Calculate the breakdown of the new curves that replace the disturbed ones.
    # The number of new curves to draw is equal to the number of disturbed curves.
    new_curves_to_draw = disturbed_curves
    
    # One-fifth of the new curves are parabolic.
    fraction_parabolic = 1/5
    new_parabolic_curves = int(new_curves_to_draw * fraction_parabolic)

    # Two-ninths are elliptical.
    fraction_elliptical = 2/9
    new_elliptical_curves = int(new_curves_to_draw * fraction_elliptical)
    
    # The rest are circular.
    new_circular_curves = new_curves_to_draw - new_parabolic_curves - new_elliptical_curves

    # Step 4: Determine the final total number of curves in the restored pattern.
    # This is the sum of the original, unchanged curves and all the new curves drawn.
    final_total_curves = unchanged_curves + new_curves_to_draw
    
    # Step 5: Print the breakdown and the final equation.
    print("Initial State:")
    print(f"Total curves in the original pattern: {total_curves}")
    print("-" * 40)
    
    print("After Disruption:")
    print(f"Number of curves that remained unchanged: {unchanged_curves}")
    print(f"Number of curves that were disturbed: {disturbed_curves}")
    print("-" * 40)

    print("Restoration - New Curves Drawn:")
    print(f"New parabolic curves drawn: {new_parabolic_curves}")
    print(f"New elliptical curves drawn: {new_elliptical_curves}")
    print(f"New circular curves drawn: {new_circular_curves}")
    print("-" * 40)
    
    print("Final Restored Pattern Calculation:")
    print("The final pattern is composed of the unchanged curves plus the new curves.")
    print("Equation: Unchanged + Parabolic + Elliptical + Circular = Final Total")
    print(f"{unchanged_curves} + {new_parabolic_curves} + {new_elliptical_curves} + {new_circular_curves} = {final_total_curves}")


solve_rangoli_puzzle()
<<<360>>>
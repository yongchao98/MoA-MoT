import math

def solve_rangoli_puzzle():
    """
    This function solves the rangoli puzzle by calculating the total number of
    curves in the final, restored pattern.
    """
    
    # Step 1: Identify the key numbers given in the problem description.
    # The original pattern had a total of 360 curves.
    initial_total_curves = 360
    
    # The poem and text specify that 90 curves remained unchanged.
    unchanged_curves = 90
    
    # Step 2: Calculate the number of curves that were altered and need to be redrawn.
    # This is the difference between the initial total and the unchanged curves.
    # This group of curves represents "those that left their plotted way."
    altered_curves = initial_total_curves - unchanged_curves
    
    # Step 3: Calculate the number of each type of new curve the master must draw.
    # The fractions for the new curve types apply to the total number of altered curves.
    # One-fifth of the altered curves become parabolic.
    parabolic_curves = math.trunc(1/5 * altered_curves)
    
    # Two-ninths of the altered curves become elliptical.
    elliptical_curves = math.trunc(2/9 * altered_curves)
    
    # The remainder of the altered curves become circular.
    circular_curves = altered_curves - (parabolic_curves + elliptical_curves)
    
    # Step 4: Calculate the total number of curves in the final, restored pattern.
    # This is the sum of the unchanged curves and all the newly drawn curves.
    # The golden ratio detail describes the nature of the curves, not their quantity.
    # Each altered curve is replaced by one new curve to ensure the pattern is complete.
    final_total_curves = unchanged_curves + parabolic_curves + elliptical_curves + circular_curves
    
    # Step 5: Print the breakdown and the final equation as requested.
    print("The final restored Rangoli pattern is composed of the following curves:")
    print(f"- Unchanged Original Curves: {unchanged_curves}")
    print(f"- New Parabolic Curves: {parabolic_curves}")
    print(f"- New Elliptical Curves: {elliptical_curves}")
    print(f"- New Circular Curves: {circular_curves}")
    print("\nThe total number of curves in the master's restored pattern is the sum of all parts:")
    print(f"{unchanged_curves} + {parabolic_curves} + {elliptical_curves} + {circular_curves} = {final_total_curves}")

solve_rangoli_puzzle()
def solve_rangoli_puzzle():
    """
    Calculates the total number of curves in the restored rangoli pattern.
    """
    # Step 1: Define the initial and unchanged number of curves based on the problem statement.
    total_initial_curves = 360
    unchanged_curves = 90

    # Step 2: Calculate the number of curves that were disturbed.
    # These are the ones that "left their plotted way".
    disturbed_curves = total_initial_curves - unchanged_curves

    # Step 3: Calculate the number of each type of new curve that replaces the disturbed ones.
    # The fractions are applied to the number of disturbed curves.
    parabolic_curves = (1/5) * disturbed_curves
    elliptical_curves = (2/9) * disturbed_curves
    
    # The rest of the disturbed curves become circular.
    circular_curves = disturbed_curves - parabolic_curves - elliptical_curves

    # Step 4: Calculate the total number of curves in the final, restored pattern.
    # This is the sum of the curves that were never disturbed and all the new curves that were drawn.
    total_final_curves = unchanged_curves + parabolic_curves + elliptical_curves + circular_curves

    # Ensure the results are integers for the final output
    unchanged_curves = int(unchanged_curves)
    parabolic_curves = int(parabolic_curves)
    elliptical_curves = int(elliptical_curves)
    circular_curves = int(circular_curves)
    total_final_curves = int(total_final_curves)

    # Step 5: Print the final equation showing the breakdown of the total curves.
    print("To restore the pattern, the master creates a new design where the total number of curves is calculated as follows:")
    print(f"Unchanged Curves: {unchanged_curves}")
    print(f"New Parabolic Curves: {parabolic_curves}")
    print(f"New Elliptical Curves: {elliptical_curves}")
    print(f"New Circular Curves: {circular_curves}")
    print("\nFinal Equation:")
    print(f"{unchanged_curves} (Unchanged) + {parabolic_curves} (Parabolic) + {elliptical_curves} (Elliptical) + {circular_curves} (Circular) = {total_final_curves}")
    
    print("\nThe total number of curves in the restored pattern is:")
    print(f"<<<{total_final_curves}>>>")

solve_rangoli_puzzle()
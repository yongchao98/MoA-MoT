def solve_rangoli_problem():
    """
    This function calculates the total number of curves the master needs to draw
    to restore the Rangoli pattern based on the provided information.
    """

    # Step 1: Establish the initial and remaining curves from the problem's explicit statements.
    initial_total_curves = 360
    undisturbed_curves = 90
    
    print(f"The original pattern had {initial_total_curves} curves.")
    print(f"The number of curves that remained undisturbed is {undisturbed_curves}.")

    # Step 2: Calculate the number of curves that were disturbed and need to be redrawn.
    disturbed_curves = initial_total_curves - undisturbed_curves
    print(f"This means the number of curves that were disturbed is {initial_total_curves} - {undisturbed_curves} = {disturbed_curves}.")
    print("-" * 20)

    # Step 3: Calculate the number of each type of new curve the master must draw.
    # These new curves replace the disturbed ones.
    print("The master must draw new curves to replace the disturbed ones.")
    
    # Parabolic curves are 1/5 of the disturbed curves.
    parabolic_curves = disturbed_curves // 5
    print(f"Number of new parabolic curves (1/5 of {disturbed_curves}): {parabolic_curves}")
    
    # Elliptical curves are 2/9 of the disturbed curves.
    elliptical_curves = (2 * disturbed_curves) // 9
    print(f"Number of new elliptical curves (2/9 of {disturbed_curves}): {elliptical_curves}")

    # The rest of the new curves are circular.
    circular_curves = disturbed_curves - (parabolic_curves + elliptical_curves)
    print(f"Number of new circular curves (the rest): {circular_curves}")
    print("-" * 20)

    # Step 4: Calculate the total number of curves the master must draw.
    # This is the sum of all the new curves.
    total_curves_to_draw = parabolic_curves + elliptical_curves + circular_curves
    
    print("The question asks for the total number of curves the master must draw.")
    print("This is the sum of all the new curves needed for the restoration.")
    print("\nFinal Equation:")
    print(f"{parabolic_curves} (parabolic) + {elliptical_curves} (elliptical) + {circular_curves} (circular) = {total_curves_to_draw}")

# Execute the function to find the answer.
solve_rangoli_problem()
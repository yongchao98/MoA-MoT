def solve_rangoli_puzzle():
    """
    Calculates the number of curves a Rangoli master must redraw based on the poem's scenario.
    """
    # Known values from the problem statement.
    undisturbed_curves = 90
    
    # Fractions of curves affected.
    # 3/8 lost shape, 1/4 found new paths.
    # We assume these are distinct groups.
    fraction_lost_shape_num = 3
    fraction_lost_shape_den = 8
    fraction_new_paths_num = 1
    fraction_new_paths_den = 4

    # Calculate the total fraction of affected curves: 3/8 + 1/4 = 5/8
    total_affected_fraction_num = 5
    total_affected_fraction_den = 8

    # The remaining fraction is undisturbed: 1 - 5/8 = 3/8
    undisturbed_fraction_num = 3
    undisturbed_fraction_den = 8

    # Calculate the total number of original curves (T).
    # (3/8) * T = 90  =>  T = 90 * (8/3)
    total_original_curves = undisturbed_curves * undisturbed_fraction_den // undisturbed_fraction_num

    # The number of curves to redraw is the total number of affected curves.
    # This can be calculated as Total Curves - Undisturbed Curves.
    curves_to_redraw = total_original_curves - undisturbed_curves

    # --- Output the reasoning and the final answer ---
    print("To solve the puzzle, we first determine the original total number of curves in the pattern.")
    print(f"The fraction of curves that 'lost shape' is {fraction_lost_shape_num}/{fraction_lost_shape_den}.")
    print(f"The fraction that 'found new paths' is {fraction_new_paths_num}/{fraction_new_paths_den}.")
    print(f"Total fraction of affected curves = {fraction_lost_shape_num}/{fraction_lost_shape_den} + {fraction_new_paths_num}/{fraction_new_paths_den} = {total_affected_fraction_num}/{total_affected_fraction_den}.")
    print(f"This means the fraction of undisturbed curves is 1 - {total_affected_fraction_num}/{total_affected_fraction_den} = {undisturbed_fraction_num}/{undisturbed_fraction_den}.")
    
    print(f"\nGiven that {undisturbed_curves} curves remained undisturbed, we can find the total:")
    print(f"Total Original Curves = {undisturbed_curves} / ({undisturbed_fraction_num}/{undisturbed_fraction_den}) = {total_original_curves}")

    print("\nThe number of curves the master must draw is the total number of affected curves.")
    print("This is the difference between the original total and the undisturbed curves.")
    
    print("\nFinal Equation:")
    print(f"{total_original_curves} (Total Curves) - {undisturbed_curves} (Undisturbed Curves) = {curves_to_redraw} (Curves to Redraw)")

solve_rangoli_puzzle()
<<<150>>>
def solve_rangoli_riddle():
    """
    Solves the Rangoli riddle by calculating the total number of curves
    in the restored pattern based on the poem's description.
    """

    # Step 1: Determine the original total number of curves.
    # The poem states 3/8 of curves "lost shape" and 1/4 "found new paths".
    # This means the fraction of curves that remained undisturbed is: 1 - 3/8 - 1/4
    frac_lost = 3/8
    frac_new_path = 1/4
    # To subtract, we find a common denominator: 1 - 3/8 - 2/8 = 3/8
    frac_remaining = 1 - frac_lost - frac_new_path

    # We are told that 90 curves remained.
    undisturbed_curves = 90

    # If (3/8) of the total is 90, we can find the original total.
    original_total_curves = int(undisturbed_curves / frac_remaining)
    print(f"First, we calculate the original total number of curves in the pattern.")
    print(f"The fraction of remaining curves is 1 - 3/8 - 1/4 = {int(frac_remaining * 8)}/8.")
    print(f"Calculation for original total: {undisturbed_curves} / ({int(frac_remaining*8)}/8) = {original_total_curves}\n")

    # Step 2: Calculate the number of curves in each group.
    lost_shape_curves = int(frac_lost * original_total_curves)
    new_path_curves = int(frac_new_path * original_total_curves)
    print("Next, we find the number of curves in each group:")
    print(f" - Curves that remained undisturbed: {undisturbed_curves}")
    print(f" - Curves that 'lost shape': 3/8 of {original_total_curves} = {lost_shape_curves}")
    print(f" - Curves that 'found new paths': 1/4 of {original_total_curves} = {new_path_curves}\n")

    # Step 3: Calculate the breakdown of the remade curves.
    # The poem breaks down "those that left their plotted way". This must be a number
    # divisible by 5 and 9. Of our groups (90, 60, 90), only 90 works.
    # So, the 90 curves that "lost shape" are the ones being remade into new types.
    remade_group_size = lost_shape_curves
    num_parabolic = int((1/5) * remade_group_size)
    num_elliptical = int((2/9) * remade_group_size)
    num_circular = remade_group_size - num_parabolic - num_elliptical
    
    print("The poem then details the new shapes for the 90 curves that 'lost shape':")
    print(f" - New Parabolic Curves: 1/5 of {remade_group_size} = {num_parabolic}")
    print(f" - New Elliptical Curves: 2/9 of {remade_group_size} = {num_elliptical}")
    print(f" - New Circular Curves (the rest): {remade_group_size} - {num_parabolic} - {num_elliptical} = {num_circular}\n")
    
    # Step 4: Calculate the total curves in the final, restored pattern.
    # The final pattern consists of the undisturbed curves plus all the replaced curves.
    total_final_curves = undisturbed_curves + lost_shape_curves + new_path_curves
    
    print("Finally, to find the total curves the master must place, we sum the parts of the restored pattern.")
    print("This includes the original undisturbed curves and all the curves that were redrawn.")
    print(f"Final Total = (Undisturbed Curves) + (Replaced 'Lost Shape' Curves) + (Replaced 'New Path' Curves)")
    print(f"Final Total = {undisturbed_curves} + {lost_shape_curves} + {new_path_curves} = {total_final_curves}")

solve_rangoli_riddle()
<<<240>>>
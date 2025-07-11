def solve_rangoli_problem():
    """
    Calculates the total number of curves the master must draw to restore the Rangoli pattern.
    """
    # Step 1: Define the known, consistent value from the problem.
    unchanged_curves = 90
    print(f"The number of curves that remained unchanged is {unchanged_curves}.")
    print("-" * 20)

    # Step 2: Determine the fraction of curves that were unchanged.
    # The total fraction of affected curves is the sum of those that lost shape (3/8)
    # and those that found new paths (1/4).
    # Fraction of affected curves = 3/8 + 1/4 = 3/8 + 2/8 = 5/8.
    # The fraction of unchanged curves is 1 - 5/8 = 3/8.
    print("The fraction of affected curves is 3/8 (lost shape) + 1/4 (found new paths) = 5/8.")
    print("Therefore, the fraction of unchanged curves is 1 - 5/8 = 3/8.")
    print("-" * 20)

    # Step 3: Calculate the original total number of curves (T).
    # We know that (3/8) * T = 90.
    # So, T = 90 * 8 / 3.
    total_original_curves = int(unchanged_curves * 8 / 3)
    print(f"The total number of original curves is calculated as: {unchanged_curves} / (3/8) = {total_original_curves}.")
    print("-" * 20)

    # Step 4: Calculate the number of curves in each affected group.
    # These are the curves the master must redraw.
    lost_shape_curves = int((3/8) * total_original_curves)
    found_new_paths_curves = int((1/4) * total_original_curves)

    # Step 5: Calculate the total number of curves to be redrawn.
    total_curves_to_draw = lost_shape_curves + found_new_paths_curves

    # Step 6: Print the final breakdown and the answer.
    print("To restore the pattern, the master must redraw all affected curves.")
    print(f"Number of curves that lost their shape: (3/8) * {total_original_curves} = {lost_shape_curves}")
    print(f"Number of curves that found new paths: (1/4) * {total_original_curves} = {found_new_paths_curves}")
    print("\nFinal Equation:")
    print(f"Total curves to draw = {lost_shape_curves} + {found_new_paths_curves} = {total_curves_to_draw}")

solve_rangoli_problem()
<<<150>>>
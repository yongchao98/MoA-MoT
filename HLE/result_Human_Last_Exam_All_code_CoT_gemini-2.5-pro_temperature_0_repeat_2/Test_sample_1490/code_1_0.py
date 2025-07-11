def solve_rangoli_puzzle():
    """
    Calculates the total number of curves in the restored Rangoli pattern
    based on the information provided in the poem and prose.
    """

    # 1. The number of curves that remained unchanged is given.
    unaffected_curves = 90

    # 2. The fractions of curves that were disturbed are given.
    fraction_lost_form = 3/8
    fraction_new_paths = 1/4

    # 3. Calculate the total fraction of curves that were affected.
    total_fraction_affected = fraction_lost_form + fraction_new_paths

    # 4. From this, calculate the fraction of curves that were unaffected.
    # This corresponds to the 90 curves that "stay true and bright".
    fraction_unaffected = 1 - total_fraction_affected

    # 5. Calculate the original total number of curves in the pattern.
    # Restoring the pattern means the final total will be the same as the original.
    total_curves = unaffected_curves / fraction_unaffected

    # 6. Calculate the number of curves in each group of replacements.
    # These curves are removed and then redrawn by the master.
    curves_lost_form_replaced = total_curves * fraction_lost_form
    curves_new_paths_replaced = total_curves * fraction_new_paths

    # 7. The final pattern is composed of the unaffected curves plus all the replaced curves.
    # The code will print the breakdown of this final composition.
    print("The final pattern is composed of three groups of curves:")
    print(f"- The curves that were unaffected: {int(unaffected_curves)}")
    print(f"- The curves that lost their form and were replaced: {int(curves_lost_form_replaced)}")
    print(f"- The curves that found new pathways and were replaced: {int(curves_new_paths_replaced)}")
    print("\nThe total number of curves in the restored pattern is the sum of these parts.")
    
    # 8. Display the final equation showing how the total is calculated.
    final_equation = (
        f"{int(unaffected_curves)} (unaffected) + "
        f"{int(curves_lost_form_replaced)} (replaced) + "
        f"{int(curves_new_paths_replaced)} (replaced) = "
        f"{int(total_curves)}"
    )
    print("\nFinal Equation:")
    print(final_equation)

    # The final answer is the total number of curves.
    print(f"\nThe total number of curves the master must draw to restore the pattern is {int(total_curves)}.")


solve_rangoli_puzzle()
<<<240>>>
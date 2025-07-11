import math

def solve_rangoli_problem():
    """
    Calculates the total number of curves in the restored rangoli pattern based on the poem's scenario.
    """
    # Step 1: Determine the original number of curves.
    # The denominators in the poem are 8, 4, 5, and 9.
    # The least common multiple (LCM) is the most logical starting total.
    # LCM(8, 4, 5, 9) = 360. This is also given in the background info.
    original_curves = 360

    # Step 2: Identify the number of curves that remain untouched.
    # The poem gives a specific condition: "If ninety curves stay true and bright".
    curves_staying_true = 90

    # Step 3: Calculate the number of new curves the master must draw.
    # The poem states "one-fourth found new paths to trace", which applies to the original total.
    new_curves_fraction_numerator = 1
    new_curves_fraction_denominator = 4
    new_curves_drawn = original_curves * new_curves_fraction_numerator / new_curves_fraction_denominator
    new_curves_drawn = int(new_curves_drawn)

    # Step 4: Calculate the total number of curves in the final, restored pattern.
    # This is the sum of the curves that remained and the new ones drawn.
    total_curves_in_restored_pattern = curves_staying_true + new_curves_drawn

    # Step 5: Output the final equation and the result.
    # The question is "how many total curves must the master draw to restore the pattern".
    # This is the final count of curves in the completed pattern.
    print("The final pattern is composed of the curves that remained plus the new curves drawn.")
    print("Final Equation:")
    print(f"{curves_staying_true} (curves that stayed true) + {new_curves_drawn} (new curves drawn) = {total_curves_in_restored_pattern}")
    print("\nThus, the total number of curves in the restored pattern is:")
    print(f"<<<{total_curves_in_restored_pattern}>>>")

solve_rangoli_problem()
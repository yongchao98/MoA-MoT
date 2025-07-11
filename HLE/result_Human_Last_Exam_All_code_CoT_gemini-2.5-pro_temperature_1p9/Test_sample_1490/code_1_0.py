import math

def solve_rangoli_puzzle():
    """
    Solves the Rangoli puzzle by calculating the total number of curves in the restored pattern.
    """
    # Step 1: Define initial total curves.
    # The problem gives two potential starting points: 360 total curves, or 90 remaining curves.
    # If 90 curves are the 5/8 that remain, the total is 144. However, calculating 1/5 of the
    # disturbed curves (54) results in a non-integer (10.8), which is not possible.
    # Therefore, we proceed with the initial total of 360, which yields whole numbers.
    total_initial_curves = 360

    # Step 2: Calculate the number of disturbed and undisturbed curves.
    # Three-eighths of the curves were disturbed.
    fraction_disturbed = 3/8
    disturbed_curves = total_initial_curves * fraction_disturbed
    undisturbed_curves = total_initial_curves - disturbed_curves

    # Step 3: Calculate the number of each new curve type drawn to replace the disturbed ones.
    # One-fifth of the disturbed curves become parabolic.
    parabolic_curves = disturbed_curves * (1/5)
    # Two-ninths of the disturbed curves become elliptical.
    elliptical_curves = disturbed_curves * (2/9)
    # The remaining disturbed curves become circular.
    circular_curves = disturbed_curves - parabolic_curves - elliptical_curves

    # Ensure all calculated numbers are integers for the final equation.
    undisturbed_curves = int(undisturbed_curves)
    parabolic_curves = int(parabolic_curves)
    elliptical_curves = int(elliptical_curves)
    circular_curves = int(circular_curves)
    
    # Step 4: Calculate the total number of curves in the restored pattern.
    # This is the sum of the undisturbed curves and all the new replacement curves.
    total_final_curves = undisturbed_curves + parabolic_curves + elliptical_curves + circular_curves

    # Step 5: Display the breakdown and the final equation.
    print("The final restored pattern is composed of several types of curves:")
    print(f"- Original Undisturbed Curves: {undisturbed_curves}")
    print(f"- New Parabolic Curves: {parabolic_curves}")
    print(f"- New Elliptical Curves: {elliptical_curves}")
    print(f"- New Circular Curves: {circular_curves}")
    print("\nThe question asks for the total number of curves in the restored pattern.")
    print("This is the sum of the original curves left untouched and all the new curves the master must draw.")
    print("\nFinal Calculation:")
    print(f"{undisturbed_curves} (Undisturbed) + {parabolic_curves} (Parabolic) + {elliptical_curves} (Elliptical) + {circular_curves} (Circular) = {total_final_curves}")

solve_rangoli_puzzle()
<<<360>>>
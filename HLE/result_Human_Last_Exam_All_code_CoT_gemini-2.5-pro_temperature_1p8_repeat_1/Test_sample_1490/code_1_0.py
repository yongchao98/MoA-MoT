import math

def solve_rangoli_puzzle():
    """
    Solves the Rangoli puzzle by calculating the number of new curves the master must draw.
    """
    # Step 1: Determine the Original Total Number of Curves.
    # 5/8 of the total curves remained unchanged, which is 90 curves.
    # (5/8) * total_original_curves = 90
    unchanged_curves = 90
    fraction_unchanged = 5/8
    total_original_curves = int(unchanged_curves / fraction_unchanged)

    # Step 2: Calculate the Number of Affected Curves.
    # These are the original curves that did not remain unchanged.
    affected_curves = total_original_curves - unchanged_curves

    # Step 3: Determine the Number of Redrawn Curves.
    # This number must be a multiple of 5 and 9, so a multiple of 45.
    # It must be a subset of the affected_curves.
    # The only multiple of 45 less than or equal to 54 is 45.
    redrawn_curves = 45

    # Step 4: Calculate the Number of Each New Curve Type.
    parabolic_curves = int((1/5) * redrawn_curves)
    elliptical_curves = int((2/9) * redrawn_curves)
    circular_curves = redrawn_curves - parabolic_curves - elliptical_curves

    # Step 5: Calculate and Display the Final Answer.
    # The question asks for the total number of NEW curves to be drawn.
    total_new_curves = redrawn_curves
    
    print(f"The original pattern had {total_original_curves} curves.")
    print(f"{unchanged_curves} curves remained, while {affected_curves} curves were affected by the disruption.")
    print(f"Of the affected curves, a total of {redrawn_curves} were redrawn into new forms.")
    print("\nThe master must place the following new curves:")
    print(f"- Parabolic Curves: {parabolic_curves}")
    print(f"- Elliptical Curves: {elliptical_curves}")
    print(f"- Circular Curves: {circular_curves}")
    print("\nThe final equation for the total number of new curves to draw is:")
    print(f"{parabolic_curves} + {elliptical_curves} + {circular_curves} = {total_new_curves}")

solve_rangoli_puzzle()
<<<45>>>
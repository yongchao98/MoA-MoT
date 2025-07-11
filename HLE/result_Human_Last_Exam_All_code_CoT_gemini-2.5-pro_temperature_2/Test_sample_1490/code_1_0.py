import fractions

def solve_rangoli_riddle():
    """
    This script calculates the total number of curves the master needs to draw
    based on the information provided in the riddle.
    """

    # 1. Define the known values from the problem.
    unaffected_curves = 90
    fraction_lost = fractions.Fraction(3, 8)
    fraction_rerouted = fractions.Fraction(1, 4)

    # 2. Calculate the total fraction of curves that were affected (lost or rerouted).
    total_affected_fraction = fraction_lost + fraction_rerouted

    # 3. The unaffected curves represent the remaining fraction of the original total.
    unaffected_fraction = 1 - total_affected_fraction

    # 4. Calculate the original total number of curves.
    # original_total = unaffected_curves / unaffected_fraction
    original_total = int(unaffected_curves / unaffected_fraction)

    # 5. Calculate the number of curves in each affected group.
    # These are the curves the master must redraw.
    curves_lost = int(original_total * fraction_lost)
    curves_rerouted = int(original_total * fraction_rerouted)

    # 6. The total number of curves to draw is the sum of the lost and rerouted curves.
    total_curves_to_draw = curves_lost + curves_rerouted

    # 7. Print the final answer, showing each number in the equation.
    print(f"The number of lost curves that need to be redrawn is: {curves_lost}")
    print(f"The number of curves that found new paths and need to be redrawn is: {curves_rerouted}")
    print("\nTo restore the pattern, the master must draw the sum of these two groups.")
    print(f"Final Equation: {curves_lost} + {curves_rerouted} = {total_curves_to_draw}")
    print("\nThe total number of curves the master must draw is:")
    print(f"<<<{total_curves_to_draw}>>>")

# Execute the function to solve the riddle.
solve_rangoli_riddle()
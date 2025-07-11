def solve_rangoli_problem():
    """
    Calculates the number of curves a master artist must draw to restore a Rangoli pattern
    based on the conditions described in the poem.
    """

    # Step 1: Identify the initial total number of curves.
    # The text states the base pattern contains 360 connection curves.
    initial_curves = 360

    # Step 2: Calculate the number of curves that were disturbed.
    # The poem says "Three-eighths of curves lost shape serene".
    fraction_disturbed = 3/8
    disturbed_curves = int(fraction_disturbed * initial_curves)

    # Step 3: Determine the number of curves the master must draw.
    # To "restore the pattern", the master must replace all the disturbed curves.
    # This brings the total count back to the original 360.
    curves_to_draw = disturbed_curves

    # Step 4: Verify this number by breaking it down into the new types.
    # The poem describes the breakdown of the newly drawn curves:
    # "One-fifth through parabolas sway"
    # "Two-ninths in ellipse arcs reside"
    # "The rest in circles now abide"
    parabolic_curves = int((1/5) * curves_to_draw)
    elliptical_curves = int((2/9) * curves_to_draw)
    circular_curves = curves_to_draw - parabolic_curves - elliptical_curves

    # The fact that these calculations result in whole numbers (1/5 * 135 = 27, 2/9 * 135 = 30)
    # confirms that 135 is the correct number of curves to be drawn.

    # Step 5: Present the final answer and the breakdown equation.
    print("The master must draw a total of {} curves to restore the pattern.".format(curves_to_draw))
    print("The breakdown of these new curves is as follows:")
    print("{} (Parabolic) + {} (Elliptical) + {} (Circular) = {} (Total Curves to Draw)".format(
        parabolic_curves, elliptical_curves, circular_curves, curves_to_draw
    ))

    # The final answer is the total number of curves to draw.
    print("\n---")
    print("Final Answer:")
    # The question asks for the total curves the master must draw.
    print("<<<{}>>>".format(curves_to_draw))

solve_rangoli_problem()
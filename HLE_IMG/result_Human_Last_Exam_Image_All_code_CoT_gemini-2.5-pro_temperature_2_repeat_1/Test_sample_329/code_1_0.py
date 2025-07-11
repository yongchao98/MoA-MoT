def solve_simulation():
    """
    Analyzes the four triple-slit simulations to find the unique parameter for each
    and constructs the corresponding four-digit code.
    """

    # Analysis for Simulation 1
    # - Brightness (h): The right slit (3) is brightest, while 1 and 2 are equal. -> h is unique for slit 3.
    # - Beam Spread (w): The spread of the beams looks identical for all three slits. -> w is identical.
    # - Fringe Spacing (u): The spacing between fringes decreases from left to right. -> u increases monotonically.
    # The unique parameter is h3.
    sim1_code = 9

    # Analysis for Simulation 2
    # - Brightness (h): The left slit (1) is brightest, while 2 and 3 are equal. -> h is unique for slit 1.
    # - Beam Spread (w): The spread of the beams looks identical for all three slits. -> w is identical.
    # - Fringe Spacing (u): The spacing between fringes increases from left to right. -> u decreases monotonically.
    # The unique parameter is h1.
    sim2_code = 7

    # Analysis for Simulation 3
    # - Brightness (h): The brightness of all three slits is identical. -> h is identical.
    # - Beam Spread (w): The beam spread decreases from left to right (widest -> narrowest), meaning slit width w increases monotonically.
    # - Fringe Spacing (u): The fringe spacing for slit 1 is larger than for slits 2 and 3, which are equal. -> u is unique for slit 1.
    # The unique parameter is u1.
    sim3_code = 4

    # Analysis for Simulation 4
    # - Brightness (h): The brightness of all three slits is identical. -> h is identical.
    # - Fringe Spacing (u): The spacing between fringes decreases from left to right. -> u increases monotonically.
    # - Beam Spread (w): The beam from slit 1 is the narrowest, while 2 and 3 are wider and equal. -> w is unique for slit 1.
    # The unique parameter is w1.
    sim4_code = 1
    
    print("Step-by-step analysis and identification of the unique parameter for each simulation:")
    print("-" * 70)
    
    print("Simulation 1:")
    print("  - Height (h): The source brightness for slit 3 is the highest. Slits 1 and 2 are equal. Unique parameter is h3.")
    print("  - The corresponding code is 9.")
    print("-" * 70)

    print("Simulation 2:")
    print("  - Height (h): The source brightness for slit 1 is the highest. Slits 2 and 3 are equal. Unique parameter is h1.")
    print("  - The corresponding code is 7.")
    print("-" * 70)

    print("Simulation 3:")
    print("  - Wave number (u): The fringe spacing for slit 1 is largest. Slits 2 and 3 are equal. Unique parameter is u1.")
    print("  - The corresponding code is 4.")
    print("-" * 70)

    print("Simulation 4:")
    print("  - Width (w): The beam spread for slit 1 is the narrowest (implying widest slit). Slits 2 and 3 are equal. Unique parameter is w1.")
    print("  - The corresponding code is 1.")
    print("-" * 70)
    
    print("Combining the codes gives the final four-digit sequence:")
    print(f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}")

solve_simulation()
<<<9741>>>
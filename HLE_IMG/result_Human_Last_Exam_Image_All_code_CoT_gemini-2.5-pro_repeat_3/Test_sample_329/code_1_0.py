def solve_simulation():
    """
    Analyzes four triple-slit simulations to find the unique parameter for each.
    
    Parameter Definitions:
    - h (height): Related to the intensity/brightness of each slit at the source.
    - w (width): Slit width, inversely proportional to the angle of the diffraction envelope (the "V" shape).
    - u (wave number): Related to the wave's phase. A linear change causes a tilt; a single different value causes asymmetry.

    Rules for each simulation:
    1. One parameter is identical across all three slits.
    2. One parameter is monotonic (strictly increasing or decreasing).
    3. One parameter has one slit that is different from the other two (the "unique" parameter).

    Encoding:
    w1->1, w2->2, w3->3
    u1->4, u2->5, u3->6
    h1->7, h2->8, h3->9
    """

    # --- Simulation 1 Analysis ---
    # h (intensity): Increases left to right (h1 < h2 < h3). This is MONOTONIC.
    # u (phase): The pattern is symmetric and not tilted. This implies u is IDENTICAL.
    # w (width): The central diffraction "V" is wider than the outer two, meaning w2 is smaller than w1 and w3, which appear equal. This is the UNIQUE parameter.
    # Unique parameter: w2. Code: 2
    sim1_code = 2

    # --- Simulation 2 Analysis ---
    # w (width): The diffraction angles for all three slits are the same. This means w is IDENTICAL.
    # u (phase): The entire pattern is tilted right, implying a linear phase gradient. This means u is MONOTONIC.
    # h (intensity): The outer slits (1 and 3) are equally bright, while the middle one (2) is dimmer. This is the UNIQUE parameter.
    # Unique parameter: h2. Code: 8
    sim2_code = 8

    # --- Simulation 3 Analysis ---
    # h (intensity): The brightness of all three slits appears equal. This means h is IDENTICAL.
    # w (width): The diffraction angle increases left-to-right (narrowest "V" on left, widest on right), so width decreases. This is MONOTONIC.
    # u (phase): The pattern is not tilted but is distorted and asymmetric in the center. This indicates the middle slit's phase is different. This is the UNIQUE parameter.
    # Unique parameter: u2. Code: 5
    sim3_code = 5

    # --- Simulation 4 Analysis ---
    # w (width): The diffraction angles for all three slits are the same. This means w is IDENTICAL.
    # h (intensity): Brightness increases from left to right. This is MONOTONIC.
    # u (phase): The pattern is asymmetric, with the main lobe shifted left. This suggests a phase difference in the rightmost slit. This is the UNIQUE parameter.
    # Unique parameter: u3. Code: 6
    sim4_code = 6

    # Print the individual results and the final combined sequence
    print(f"The unique parameter code for simulation 1 is: {sim1_code}")
    print(f"The unique parameter code for simulation 2 is: {sim2_code}")
    print(f"The unique parameter code for simulation 3 is: {sim3_code}")
    print(f"The unique parameter code for simulation 4 is: {sim4_code}")
    print(f"Final four-digit sequence: {sim1_code}{sim2_code}{sim3_code}{sim4_code}")

solve_simulation()
<<<2856>>>
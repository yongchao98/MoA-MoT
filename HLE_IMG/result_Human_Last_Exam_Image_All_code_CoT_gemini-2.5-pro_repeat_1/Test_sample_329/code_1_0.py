def solve_simulation():
    """
    Analyzes the four triple-slit simulations to find the unique parameter for each.
    """
    # Encoding map for parameters: w(1-3), u(4-6), h(7-9) for slits 1, 2, 3
    # Not used in code logic, but for reference.
    encoding = {
        'w1': 1, 'w2': 2, 'w3': 3,
        'u1': 4, 'u2': 5, 'u3': 6,
        'h1': 7, 'h2': 8, 'h3': 9
    }

    print("Step-by-step analysis of the simulations:")
    print("-" * 40)

    # --- Simulation 1 Analysis ---
    print("Analysis of Simulation 1:")
    print("1. Brightness (h): The intensity of the sources at the bottom increases from left to right (h1 < h2 < h3). This follows Rule 2 (monotonic).")
    print("2. Beam Spread (w): The angular spread of the three beams appears to be the same. This follows Rule 1 (identical).")
    print("3. By elimination, Wave Number (u) must follow Rule 3 (unique). The entire pattern is bent to the left. This asymmetry indicates a phase difference, and a bend to the left suggests the rightmost slit (slit 3) is the unique one.")
    print("Conclusion: The unique parameter is u3.")
    sim1_code = 6
    print(f"The code for Simulation 1 is {sim1_code}.\n")

    # --- Simulation 2 Analysis ---
    print("Analysis of Simulation 2:")
    print("1. Brightness (h): The intensity of the three sources appears to be identical. This follows Rule 1 (identical).")
    print("2. Beam Spread (w): The spread of the beams decreases from left to right (leftmost is widest, rightmost is narrowest). This implies w1 > w2 > w3. This follows Rule 2 (monotonic).")
    print("3. By elimination, Wave Number (u) must follow Rule 3 (unique). The pattern is bent to the right, suggesting the leftmost slit (slit 1) is the unique one.")
    print("Conclusion: The unique parameter is u1.")
    sim2_code = 4
    print(f"The code for Simulation 2 is {sim2_code}.\n")

    # --- Simulation 3 Analysis ---
    print("Analysis of Simulation 3:")
    print("1. Brightness (h): The outer two slits are bright and have similar intensity, while the middle slit is significantly dimmer. This means h1 = h3, and h2 is different. This follows Rule 3 (unique).")
    print("Conclusion: The unique parameter is h2.")
    sim3_code = 8
    print(f"The code for Simulation 3 is {sim3_code}.\n")

    # --- Simulation 4 Analysis ---
    print("Analysis of Simulation 4:")
    print("1. Pattern Symmetry (u): The overall interference pattern is symmetric with no left or right bending. This indicates the wave number u is the same for all three slits. This follows Rule 1 (identical).")
    print("2. Brightness (h): The intensity of the sources increases from left to right (h1 < h2 < h3). This follows Rule 2 (monotonic).")
    print("3. By elimination, Slit Width (w) must follow Rule 3 (unique). The beam from the middle slit is visibly more spread out than the beams from the outer two slits, which look similar. A wider spread means a narrower slit. Thus, w2 is unique (and smaller than w1 and w3).")
    print("Conclusion: The unique parameter is w2.")
    sim4_code = 2
    print(f"The code for Simulation 4 is {sim4_code}.\n")

    # --- Final Result ---
    print("-" * 40)
    print("The final four-digit sequence is formed by concatenating the codes from each simulation.")
    print(f"The sequence is: {sim1_code}{sim2_code}{sim3_code}{sim4_code}")

solve_simulation()
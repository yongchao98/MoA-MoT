# This script determines the four-digit sequence by analyzing the four triple-slit simulations.

def solve_simulation():
    """
    Analyzes the four simulations to find the unique parameter for each,
    then constructs the final four-digit code.
    """
    
    # Encoding map for parameters:
    # w1->1, w2->2, w3->3, u1->4, u2->5, u3->6, h1->7, h2->8, h3->9
    encoding = {
        'w1': 1, 'w2': 2, 'w3': 3,
        'u1': 4, 'u2': 5, 'u3': 6,
        'h1': 7, 'h2': 8, 'h3': 9
    }

    # --- Simulation 1 Analysis ---
    # Rule 1 (Identical): Beam tilt 'u' is the same for all (symmetric pattern).
    # Rule 2 (Monotonic): Slit brightness 'h' increases left to right.
    # Rule 3 (Unique): Central beam width is different from the outer two, so 'w2' is unique.
    sim1_unique_param = 'w2'
    sim1_code = encoding[sim1_unique_param]

    # --- Simulation 2 Analysis ---
    # Rule 1 (Identical): Beam widths 'w' are the same for all slits.
    # Rule 2 (Monotonic): The entire pattern is shifted, indicating 'u' changes monotonically.
    # Rule 3 (Unique): The right slit 'h3' has a different brightness from the other two.
    sim2_unique_param = 'h3'
    sim2_code = encoding[sim2_unique_param]

    # --- Simulation 3 Analysis ---
    # Rule 1 (Identical): Slit brightness 'h' is the same for all slits.
    # Rule 2 (Monotonic): Beam widths 'w' are all different, indicating 'w' changes monotonically.
    # Rule 3 (Unique): The central beam tilt 'u2' is different from the outer two.
    sim3_unique_param = 'u2'
    sim3_code = encoding[sim3_unique_param]

    # --- Simulation 4 Analysis ---
    # Rule 1 (Identical): Beam widths 'w' are the same for all slits.
    # Rule 2 (Monotonic): Slit brightness 'h' decreases from left to right.
    # Rule 3 (Unique): The right beam tilt 'u3' is different from the other two.
    sim4_unique_param = 'u3'
    sim4_code = encoding[sim4_unique_param]

    # Combine the codes into the final sequence
    final_code = f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}"

    print(f"Code for Simulation 1 (Unique parameter {sim1_unique_param}): {sim1_code}")
    print(f"Code for Simulation 2 (Unique parameter {sim2_unique_param}): {sim2_code}")
    print(f"Code for Simulation 3 (Unique parameter {sim3_unique_param}): {sim3_code}")
    print(f"Code for Simulation 4 (Unique parameter {sim4_unique_param}): {sim4_code}")
    print(f"\nFinal four-digit sequence: {final_code}")

solve_simulation()
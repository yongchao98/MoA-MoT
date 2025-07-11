def solve_slit_simulations():
    """
    This script details the analysis of four triple-slit simulations
    to determine the unique parameter for each and construct the final code.
    """

    # Parameter codes:
    # w1:1, w2:2, w3:3
    # u1:4, u2:5, u3:6
    # h1:7, h2:8, h3:9

    # Analysis Results
    # Sim 1: h is monotonic (intensity L->R increase), w is identical (beam widths equal).
    # -> u is unique. Pattern is distorted in the center, so u2 is unique.
    sim1_code = 5

    # Sim 2: h is identical (intensities equal), w is monotonic (beam widths L->R decrease).
    # -> u is unique. Pattern is pulled left, so u1 is unique.
    sim2_code = 4

    # Sim 3: u is monotonic (pattern tilts right), w is identical (beam widths equal).
    # -> h is unique. Middle source is dimmer, so h2 is unique.
    sim3_code = 8

    # Sim 4: h is identical (intensities equal), u is monotonic (pattern tilts left).
    # -> w is unique. Right beam is wider, so w3 is unique.
    sim4_code = 3

    # Constructing the final answer as an "equation" of its components
    print("The unique parameter codes for each simulation are:")
    print(f"Simulation 1: {sim1_code}")
    print(f"Simulation 2: {sim2_code}")
    print(f"Simulation 3: {sim3_code}")
    print(f"Simulation 4: {sim4_code}")
    print("-" * 20)
    
    final_sequence = f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}"
    print("The final equation is the concatenation of these codes:")
    print(f"{sim1_code} | {sim2_code} | {sim3_code} | {sim4_code}  =>  {final_sequence}")

solve_slit_simulations()
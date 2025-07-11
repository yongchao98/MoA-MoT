def solve_roller_puzzle():
    """
    Analyzes and solves the roller drive matching puzzle by systematically pairing
    roller configurations with their displacement curves based on geometric features.
    """
    print("### Roller Drive Puzzle Analysis ###\n")
    print("Matching strategy:")
    print("1. Match by the number of lobes on the driver (green roller), which equals the number of 'wiggles' in the plot.")
    print("2. For configurations with the same number of lobes, differentiate based on the expected speed variation (wiggle size).\n")

    # Step-by-step matching logic
    print("--- Matching Process ---\n")

    # Match based on unique lobe counts first
    # A driver with N lobes corresponds to a plot with N wiggles.
    # Lobe counts for drivers 1-8: {1:4, 2:5, 3:2, 4:4, 5:6, 6:4, 7:6, 8:3}
    # Wiggle counts for plots A-H: {A:4, B:5, C:2, D:4, E:6, F:6, G:3, H:4}

    # Plot C has 2 wiggles -> Matches the only 2-lobed driver (Config 3)
    match_C = 3
    print(f"Plot C has 2 wiggles per cycle, matching Configuration {match_C} (2-lobed driver).")

    # Plot G has 3 wiggles -> Matches the only 3-lobed driver (Config 8)
    match_G = 8
    print(f"Plot G has 3 wiggles per cycle, matching Configuration {match_G} (3-lobed driver).")

    # Plot B has 5 wiggles -> Matches the only 5-lobed driver (Config 2)
    match_B = 2
    print(f"Plot B has 5 wiggles per cycle, matching Configuration {match_B} (5-lobed driver).")

    # Analyze the 6-lobe configurations (5, 7) and 6-wiggle plots (E, F)
    print("\n--- Group: 6 Lobes/Wiggles ---")
    # Config 7 has much more extreme/pointy shapes than Config 5, implying greater speed variation.
    # Plot E has much larger wiggles (greater variation) than Plot F.
    match_E = 7
    match_F = 5
    print(f"Config 7 is more extreme than 5. Plot E has larger wiggles than F.")
    print(f"Therefore, Plot E = {match_E} and Plot F = {match_F}.")

    # Analyze the 4-lobe configurations (1, 4, 6) and 4-wiggle plots (A, D, H)
    print("\n--- Group: 4 Lobes/Wiggles ---")
    # Rank by smoothness/expected variation: Config 6 (smoothest) < Config 1 < Config 4 (most extreme)
    # Rank by wiggle size: Plot H (smallest) < Plot A < Plot D
    match_H = 6
    print(f"Config 6 is the smoothest and will produce the least variation.")
    print(f"Plot H has the smallest wiggles. Therefore, Plot H = {match_H}.")
    
    # Differentiate between 1 and 4. Config 4 has a 4-lobed driver and a 2-lobed driven.
    # This 4:2 symmetry results in a repeating pattern of wiggles (paired wiggles).
    # Plot D shows this paired-wiggle pattern. Plot A's wiggles are more uniform.
    match_D = 4
    match_A = 1
    print(f"Config 4's 4:2 symmetry causes the paired wiggles seen in Plot D. Therefore, Plot D = {match_D}.")
    print(f"By elimination, the remaining 4-wiggle plot A matches Config {match_A}.")
    
    # Final Summary
    print("\n### Final Pairings ###")
    pairings = {
        'A': match_A, 'B': match_B, 'C': match_C, 'D': match_D,
        'E': match_E, 'F': match_F, 'G': match_G, 'H': match_H
    }
    for plot, config in sorted(pairings.items()):
        print(f"Plot {plot} = Configuration {config}")
        
    final_sequence = "".join(str(pairings[key]) for key in sorted(pairings.keys()))
    
    print("\n-------------------------------------------------------------")
    print("The final sequence for plots A-H is:")
    print(final_sequence)
    print("-------------------------------------------------------------")

solve_roller_puzzle()
def solve_tapestry_puzzle():
    """
    This script systematically analyzes the possible scenarios to identify
    the mutated enzyme and the original color of the orange patch.
    """
    print("Analyzing the metabolic pathways and mutations...\n")

    # Scenario 1: The original patch was yellow.
    print("--- Testing Hypothesis: Original Colour was YELLOW ---")

    # Lacking Enzyme A
    print("If Enzyme A is missing:")
    print("  The first step (Yellow -> Red) is blocked.")
    print("  Result: The patch remains YELLOW.\n")

    # Lacking Enzyme B
    print("If Enzyme B is missing:")
    print("  Enzyme A converts the yellow pigment to a red intermediate.")
    print("  The pathway stops there, and the red intermediate accumulates.")
    print("  The mixture of the remaining yellow pigment and the new red intermediate appears ORANGE.")
    print("  Result: This matches the observation!\n")

    # Lacking Enzyme C
    print("If Enzyme C is missing:")
    print("  The pathway proceeds to the blue intermediate (Yellow -> Red -> Blue).")
    print("  The pathway stops there, and the blue intermediate accumulates.")
    print("  Result: The patch turns BLUE.\n")

    # Lacking Enzyme D
    print("If Enzyme D is missing:")
    print("  Enzyme D is not involved in the yellow pigment pathway.")
    print("  The pathway completes to the final colourless product.")
    print("  Result: The patch becomes COLOURLESS.\n")

    # Scenario 2: The original patch was blue.
    print("--- Testing Hypothesis: Original Colour was BLUE ---")

    # Lacking Enzyme A, B, or C
    print("If Enzyme A, B, or C is missing:")
    print("  These enzymes are not involved in the blue pigment pathway.")
    print("  The pathway (Blue -> Yellow) proceeds normally using Enzyme D.")
    print("  Result: The patch turns YELLOW.\n")

    # Lacking Enzyme D
    print("If Enzyme D is missing:")
    print("  The degradation of blue pigment (Blue -> Yellow) is blocked.")
    print("  Result: The patch remains BLUE.\n")

    # Conclusion
    print("--- Conclusion ---")
    print("The only scenario that results in an orange patch is when the original colour was yellow and the microbe was lacking Enzyme B.")
    
    enzyme = "B"
    original_colour = "yellow"
    
    print("\nFinal Answer:")
    # The final answer format is <enzyme>-<colour>
    print(f"{enzyme}-{original_colour}")

solve_tapestry_puzzle()
<<<B-yellow>>>
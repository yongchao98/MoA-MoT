def analyze_raphidoptera_diet():
    """
    Analyzes the dietary habits of adult Raphidiopterans to answer the user's question.
    """
    # Step 1: Define known dietary facts for adult Raphidiopterans.
    primary_diet = "Predators of small, soft-bodied arthropods (e.g., aphids, psyllids, small caterpillars, insect eggs)."
    supplementary_diet = "Pollen and nectar."

    print("Analyzing the provided choices for the diet of adult Raphidiopterans:\n")

    # Step 2: Evaluate each choice.
    choices = {
        'A': "Nectar",
        'B': "Māhoe pollen",
        'C': "Fungus",
        'D': "Karamū leaf tissue",
        'E': "Totara Aphids",
        'F': "A and E",
        'G': "D and E"
    }

    print(f"Choice A: {choices['A']}")
    print(f"Analysis: Correct. Adult snakeflies are known to supplement their diet with nectar. It provides a source of energy.")
    print("-" * 20)

    print(f"Choice E: {choices['E']}")
    print(f"Analysis: Correct. Adult snakeflies are primarily predators, and aphids are a common food source. The type of plant the aphid is on (Totara) is irrelevant to the snakefly.")
    print("-" * 20)
    
    print(f"Choice B, C, D:")
    print(f"Analysis: Incorrect. Snakeflies are not found in New Zealand (where Māhoe, Karamū, and Totara are native). More importantly, they do not feed on fungus or leaf tissue.")
    print("-" * 20)

    # Step 3: Determine the best fit.
    print("Conclusion:")
    print("Since both A (Nectar) and E (Aphids) have been recorded as food sources for adult Raphidiopterans, the most complete and accurate answer is the one that includes both.")
    print("\nThe final choice is F, which combines A and E.")

analyze_raphidoptera_diet()
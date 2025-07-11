def analyze_snakefly_diet():
    """
    Analyzes the provided options to determine the diet of adult Raphidiopterans.
    """
    # Step 1: Define known dietary facts for adult Raphidiopterans (Snakeflies).
    # They are predators and also feed on pollen/nectar.
    primary_diet = "small, soft-bodied arthropods (e.g., aphids, mites)"
    supplementary_diet = "pollen and nectar"

    print("Analyzing the diet of adult Raphidiopterans...\n")

    # Step 2: Evaluate each option.
    # Option A: Nectar
    analysis_A = "Nectar" in supplementary_diet
    print(f"1. Is Nectar part of their diet? {'Yes' if analysis_A else 'No'}.")
    print(f"   - Fact: Raphidiopterans are known to feed on {supplementary_diet}.\n")

    # Option B: Māhoe pollen
    # While they eat pollen, Māhoe is a New Zealand plant, and snakeflies are not native there.
    # This is likely a distractor, but the food type "pollen" is correct.
    analysis_B = "pollen" in supplementary_diet
    print(f"2. Is Māhoe pollen part of their diet? {'Partially' if analysis_B else 'No'}.")
    print(f"   - Fact: While they eat pollen, the specific mention of 'Māhoe' is a distractor, as Raphidiopterans are not native to New Zealand.\n")

    # Option C: Fungus
    analysis_C = False
    print(f"3. Is Fungus part of their diet? {'No' if not analysis_C else 'Yes'}.")
    print("   - Fact: They are not known to be fungivores.\n")

    # Option D: Karamū leaf tissue
    analysis_D = False
    print(f"4. Is Karamū leaf tissue part of their diet? {'No' if not analysis_D else 'Yes'}.")
    print("   - Fact: They are predators and pollen/nectar feeders, not herbivores that eat leaf tissue.\n")

    # Option E: Totara Aphids
    analysis_E = "aphids" in primary_diet
    print(f"5. Are Totara Aphids part of their diet? {'Yes' if analysis_E else 'No'}.")
    print(f"   - Fact: Raphidiopterans are predators that feed on {primary_diet}. Aphids are a key food source.\n")

    # Step 3: Conclude based on the analysis.
    print("-----------------------------------------")
    print("Conclusion:")
    if analysis_A and analysis_E:
        print("Both Nectar (A) and Aphids (E) are recorded food sources for adult Raphidiopterans.")
        print("Therefore, the most complete answer is F, which includes both A and E.")
        final_answer = "F"
    elif analysis_D and analysis_E:
        final_answer = "G"
    else:
        # Fallback for single answers if no combination matches.
        final_answer = "Inconclusive based on combinations"

    print(f"\nThe correct option is {final_answer}.")

# Run the analysis
analyze_snakefly_diet()
<<<F>>>
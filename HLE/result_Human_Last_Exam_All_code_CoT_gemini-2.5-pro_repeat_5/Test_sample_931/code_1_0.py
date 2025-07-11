def solve_raphidioptera_diet():
    """
    Analyzes the diet of adult Raphidiopterans (snakeflies) to determine the correct answer.
    """

    # Step 1: Define known dietary facts about adult Raphidiopterans.
    # Fact 1: They are primarily predators of small, soft-bodied insects.
    # Fact 2: They are known to supplement their diet with nectar and pollen.
    # Fact 3: Raphidioptera are not native to New Zealand.
    fact_predatory = "Adult snakeflies prey on soft-bodied insects like aphids."
    fact_supplementary = "Adult snakeflies also feed on nectar and pollen."
    fact_geographic = "Snakeflies are not found in New Zealand."

    # Step 2: Evaluate the answer choices based on these facts.
    analysis = {
        'A. Nectar': 'Correct. This is a known supplementary food source.',
        'B. Māhoe pollen': 'Incorrect. Māhoe is a New Zealand plant, and snakeflies are not found there.',
        'C. Fungus': 'Incorrect. They are not known to be fungivores.',
        'D. Karamū leaf tissue': 'Incorrect. Karamū is a New Zealand plant, and snakeflies are not herbivores.',
        'E. Totara Aphids': 'Partially correct. They eat aphids, but the mention of Totara (a NZ tree) is misleading. However, "Aphids" as a food type is correct.',
    }

    # Step 3: Conclude based on the analysis.
    # Both Nectar (A) and Aphids (a key part of E) are correct food sources.
    # Therefore, the option combining both is the most accurate.
    conclusion = "Based on the analysis, adult Raphidiopterans feed on both nectar and aphids. Therefore, option F, which combines A and E, is the most complete and correct answer."

    print("--- Analysis of Adult Raphidiopteran Diet ---")
    print(f"Fact 1: {fact_predatory}")
    print(f"Fact 2: {fact_supplementary}")
    print(f"Fact 3: {fact_geographic}")
    print("\n--- Evaluation of Choices ---")
    for choice, reason in analysis.items():
        print(f"- {choice}: {reason}")
    print("\n--- Conclusion ---")
    print(conclusion)
    
    final_answer_choice = 'F'
    print(f"\nFinal Answer Choice: {final_answer_choice}")

solve_raphidioptera_diet()
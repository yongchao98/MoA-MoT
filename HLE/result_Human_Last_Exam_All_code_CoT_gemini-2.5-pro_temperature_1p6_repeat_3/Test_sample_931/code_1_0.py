import collections

def solve_raphidioptera_diet():
    """
    Analyzes the dietary options for adult Raphidiopterans based on known facts.
    """
    # Known facts about adult Raphidiopteran diet
    known_diet_facts = {
        "predatory": True,
        "prey_type": "small, soft-bodied arthropods (e.g., aphids, psyllids)",
        "supplemental_diet": ["nectar", "pollen"]
    }

    # The provided answer choices
    choices = {
        'A': 'Nectar',
        'B': 'Māhoe pollen',
        'C': 'Fungus',
        'D': 'Karamū leaf tissue',
        'E': 'Totara Aphids',
        'F': ['A', 'E'],
        'G': ['D', 'E']
    }

    print("Analyzing the diet of adult Raphidiopterans...\n")

    # Evaluate each option
    is_correct = collections.defaultdict(bool)

    # A: Nectar
    if 'nectar' in known_diet_facts["supplemental_diet"]:
        print("Evaluating A (Nectar): Correct. Adult snakeflies are known to feed on nectar.")
        is_correct['A'] = True
    else:
        print("Evaluating A (Nectar): Incorrect.")

    # B: Māhoe pollen
    if 'pollen' in known_diet_facts["supplemental_diet"]:
        # While specific plant pollen might vary, pollen in general is correct.
        print("Evaluating B (Māhoe pollen): Plausible. Adult snakeflies feed on pollen.")
        is_correct['B'] = True
    else:
        print("Evaluating B (Māhoe pollen): Incorrect.")

    # C: Fungus
    print("Evaluating C (Fungus): Incorrect. Snakeflies are not known to be fungivores.")
    is_correct['C'] = False

    # D: Karamū leaf tissue
    print("Evaluating D (Karamū leaf tissue): Incorrect. Snakeflies are not herbivores and do not eat leaf tissue.")
    is_correct['D'] = False

    # E: Totara Aphids
    if known_diet_facts["predatory"] and "aphids" in known_diet_facts["prey_type"]:
        print("Evaluating E (Totara Aphids): Correct. Adult snakeflies are predators that feed on aphids.")
        is_correct['E'] = True
    else:
        print("Evaluating E (Totara Aphids): Incorrect.")

    print("\n---Conclusion---")
    # Evaluate combined options
    # F combines A and E.
    if is_correct['A'] and is_correct['E']:
        print("Option F combines A and E. Both Nectar (A) and Aphids (E) are recorded food sources.")
        print("This makes F the most comprehensive and correct answer.")
        final_answer = 'F'
    # G combines D and E.
    elif is_correct['D'] and is_correct['E']:
        final_answer = 'G'
    else:
        # Fallback if logic for combined choices fails
        print("Reviewing individual answers to find the best fit.")
        final_answer = None # Should not be reached in this case

    print(f"\nFinal analysis points to choice {final_answer} as the correct answer.")

solve_raphidioptera_diet()
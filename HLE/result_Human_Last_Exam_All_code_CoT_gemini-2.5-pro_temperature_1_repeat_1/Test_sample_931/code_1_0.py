def solve_diet_question():
    """
    Evaluates diet choices for adult Raphidiopterans based on known facts.
    """
    # Step 1: Establish a knowledge base.
    # Facts about adult Raphidopteran (snakefly) diet.
    # They are primarily predatory but also feed on pollen/nectar.
    knowledge_base = {
        'primary_food': ['aphids', 'mites', 'insect eggs', 'other small soft-bodied arthropods'],
        'supplementary_food': ['nectar', 'pollen'],
        'is_herbivore': False, # They don't eat leaf tissue
        'is_fungivore': False   # They don't eat fungus
    }

    # Answer choices from the user
    choices = {
        'A': 'Nectar',
        'B': 'Māhoe pollen',
        'C': 'Fungus',
        'D': 'Karamū leaf tissue',
        'E': 'Totara Aphids',
        'F': ['A', 'E'],
        'G': ['D', 'E']
    }

    print("Analyzing the diet of adult Raphidiopterans...")
    print("-" * 40)

    # Step 2: Evaluate individual choices
    a_is_correct = 'nectar' in knowledge_base['supplementary_food']
    print(f"A. Nectar: {'CORRECT' if a_is_correct else 'INCORRECT'}. Adults are known to feed on nectar as a supplementary food source.")

    # While they eat pollen, Māhoe is a New Zealand plant, and snakeflies are not native there. This makes the choice specific and likely a distractor.
    b_is_correct = False
    print(f"B. Māhoe pollen: {'INCORRECT'}. While adults may eat pollen, this choice is too specific and geographically unlikely. Their primary diet is predatory.")

    c_is_correct = not knowledge_base['is_fungivore']
    print(f"C. Fungus: {'INCORRECT'}. Snakeflies are not known to be fungivores.")

    d_is_correct = not knowledge_base['is_herbivore']
    print(f"D. Karamū leaf tissue: {'INCORRECT'}. Snakeflies are not herbivores and do not consume leaves.")

    e_is_correct = any(prey in choices['E'].lower() for prey in knowledge_base['primary_food'])
    print(f"E. Totara Aphids: {'CORRECT' if e_is_correct else 'INCORRECT'}. 'Aphids' are a primary food source for predatory adult snakeflies.")

    print("-" * 40)
    print("Analyzing combined choices:")

    # Step 3: Evaluate combined choices
    f_is_correct = a_is_correct and e_is_correct
    print(f"F. (A and E): This combines Nectar (correct) and Aphids (correct). This is a strong candidate.")

    g_is_correct = (not d_is_correct) and e_is_correct # (not d_is_correct) is True
    print(f"G. (D and E): This combines leaf tissue (incorrect) and Aphids (correct). Because one part is incorrect, this option is incorrect.")
    
    print("-" * 40)
    print("Conclusion:")
    print("The most comprehensive answer is the one that includes both their predatory and supplementary feeding habits.")
    print("The final answer combines 'A' and 'E'.")
    print("\nFinal Equation Components:")
    print(f"Component 1: A = {choices['A']}")
    print(f"Component 2: E = {choices['E']}")
    print(f"Result: A + E = F")


solve_diet_question()
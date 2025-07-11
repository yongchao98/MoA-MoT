def find_snakefly_diet():
    """
    This script determines the correct diet for adult Raphidiopterans (snakeflies)
    from a list of choices and prints the reasoning.
    """
    
    # Define the choices provided in the question.
    choices = {
        'A': 'Nectar',
        'E': 'Totara Aphids',
        'F': 'A and E'
    }

    # Step 1: Analyze the known diet of adult Raphidiopterans.
    analysis = [
        "Adult Raphidiopterans (snakeflies) are primarily predators.",
        "They are well-documented to feed on small, soft-bodied arthropods, including aphids.",
        f"This makes choice E, '{choices['E']}', a correct statement.",
        "\nFurthermore, some species of adult snakeflies are also known to supplement their diet.",
        "They have been recorded feeding on pollen and nectar from flowers.",
        f"This makes choice A, '{choices['A']}', also a correct statement."
    ]

    print("Step-by-step analysis of Raphidiopteran diet:")
    for step in analysis:
        print(f"- {step}")

    # Step 2: Conclude the best answer based on the analysis.
    conclusion = f"\nSince both statements A and E are correct, the best and most complete answer is the one that combines them."
    
    print(conclusion)
    print(f"Final Answer Choice: {choices['F']}")

find_snakefly_diet()
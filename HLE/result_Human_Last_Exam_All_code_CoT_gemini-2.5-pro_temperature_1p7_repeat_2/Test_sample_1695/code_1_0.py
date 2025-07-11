def solve_sentence_puzzle():
    """
    Solves the puzzle by generating all possible sentences and filtering for those
    that are true and meaningful.
    """
    subjects = ['Satyendranath Bose', 'The conversation', 'Sushi', 'Our Tuna']
    verbs = ['baptized', 'turned', 'contains', 'is a']
    objects = ['Bosons', 'playful', 'fish', 'Bengalese']

    # This dictionary will hold the logic for determining if a sentence is true and meaningful.
    # The keys are verbs, and the values are functions that take a subject and object
    # and return True if the resulting sentence is valid.
    # Note on Premises:
    # 1. 'Satyendranath Bose baptized Bosons' -> Physics context is relevant. 'baptized' is metaphorical for 'named'.
    # 2. 'The conversation turned playful' -> 'turned' means 'became'.
    # 3. 'Sushi contains fish' -> 'contains' implies composition.
    # 4. 'Our Tuna is a Bengalese' -> Establishes 'Our Tuna' is a Bengalese finch, not a fish.

    validation_rules = {
        'baptized': lambda s, o: s == 'Satyendranath Bose' and o == 'Bosons',

        'turned': lambda s, o: o == 'playful' and s in ['Satyendranath Bose', 'The conversation', 'Our Tuna'],

        'contains': lambda s, o:
            # Physics truth: Any physical object contains bosons (e.g., photons).
            # The conversation is abstract and does not contain them.
            (o == 'Bosons' and s in ['Satyendranath Bose', 'Sushi', 'Our Tuna']) or
            # Common knowledge truth.
            (s == 'Sushi' and o == 'fish'),

        'is a': lambda s, o:
            # Grammatically, the object must be a singular noun. 'Bosons' is plural, 'playful' is an adjective.
            # 'Satyendranath Bose is a Bengalese' is true (he was Bengali).
            (s == 'Satyendranath Bose' and o == 'Bengalese') or
            # 'Our Tuna is a Bengalese' is true by the premise.
            (s == 'Our Tuna' and o == 'Bengalese')
    }

    true_sentences = []
    counts_per_verb = {v: 0 for v in verbs}

    # Iterate through all combinations of subjects, verbs, and objects
    for verb in verbs:
        for subject in subjects:
            for obj in objects:
                # Apply the validation rule for the current verb
                if validation_rules[verb](subject, obj):
                    sentence = f"{subject} {verb} {obj}."
                    true_sentences.append(sentence)
                    counts_per_verb[verb] += 1
    
    # Print the results
    print("Found the following true and meaningful sentences:")
    for sentence in sorted(true_sentences):
        print(f"- {sentence}")
    
    print("\nCalculation:")
    count_strings = [str(c) for c in counts_per_verb.values()]
    total_count = sum(counts_per_verb.values())
    
    # The final equation as requested
    print(f"{' + '.join(count_strings)} = {total_count}")

solve_sentence_puzzle()
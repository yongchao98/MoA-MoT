def solve_sentence_puzzle():
    """
    Solves the sentence swapping puzzle by identifying all true and meaningful sentences
    that can be formed by permuting subjects and objects.
    """
    subjects = ['Satyendranath Bose', 'The conversation', 'Sushi', 'Our Tuna']
    verbs_and_predicates = {
        'baptized': 'baptized',
        'turned': 'turned',
        'contains': 'contains',
        'is': 'is'
    }
    objects = ['Bosons', 'playful', 'fish', 'a Bengalese']

    # Based on the reasoning, these are the combinations that form
    # true and meaningful sentences.
    true_combinations = {
        ('Satyendranath Bose', 'baptized', 'Bosons'),
        ('The conversation', 'turned', 'playful'),
        ('Sushi', 'contains', 'fish'),
        ('Our Tuna', 'is', 'a Bengalese'),
        ('Our Tuna', 'turned', 'playful'),
        ('Our Tuna', 'contains', 'fish'),
        ('Our Tuna', 'is', 'playful'),
        ('Satyendranath Bose', 'is', 'a Bengalese'),
        ('The conversation', 'is', 'playful'),
    }

    count = 0
    equation_parts = []
    
    print("Evaluating all possible sentences...\n")
    print("Found the following true and meaningful sentences:")
    
    # Iterate through all combinations
    for verb_key in verbs_and_predicates:
        for subj in subjects:
            for obj in objects:
                # Check if this combination is in our set of true sentences
                if (subj, verb_key, obj) in true_combinations:
                    count += 1
                    equation_parts.append("1")
                    
                    # Handle grammar for "is a" vs "is"
                    verb_phrase = verbs_and_predicates[verb_key]
                    if verb_key == 'is' and obj != 'playful':
                        sentence = f"{count}. {subj} {verb_phrase} {obj}."
                    elif verb_key == 'is' and obj == 'playful':
                         sentence = f"{count}. {subj} {verb_phrase} {obj}."
                    else:
                        sentence = f"{count}. {subj} {verb_phrase} {obj}."
                    
                    # Correction for grammar with "is a Bengalese"
                    if sentence.endswith("is a Bengalese."):
                        pass # correct
                    elif sentence.endswith(" is a fish.") or sentence.endswith(" is a Bosons."):
                        sentence = sentence.replace(" is a ", " is ") # Better grammar

                    print(sentence)

    # Print the final equation and the total count
    final_equation = " + ".join(equation_parts)
    print("\nThe total number of true and meaningful sentences is:")
    print(f"{final_equation} = {count}")


solve_sentence_puzzle()
<<<9>>>
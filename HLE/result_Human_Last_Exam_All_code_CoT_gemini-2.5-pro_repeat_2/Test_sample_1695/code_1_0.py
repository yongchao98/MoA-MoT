def solve_sentence_puzzle():
    """
    Solves the sentence puzzle by generating all possible sentences
    and filtering for those that are true and meaningful.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # These are the (subject_idx, verb_idx, object_idx) tuples
    # that form true and meaningful sentences based on common knowledge.
    true_combinations = [
        # Original sentences
        (0, 0, 0),  # Satyendranath Bose baptized Bosons.
        (1, 1, 1),  # The conversation turned playful.
        (2, 2, 2),  # Sushi contains fish.
        (3, 3, 3),  # Our Tuna is a Bengalese.

        # New true and meaningful sentences
        # --- Verb: turned ---
        (3, 1, 1),  # Our Tuna turned playful. (A cat can become playful)
        
        # --- Verb: contains ---
        (1, 2, 0),  # The conversation contains Bosons. (A conversation can be about a topic)
        (1, 2, 2),  # The conversation contains fish. (A conversation can be about a topic)
        (1, 2, 3),  # The conversation contains a Bengalese. (A conversation can be about a topic)
        (3, 2, 2),  # Our Tuna contains fish. (Cat food often contains fish)
        
        # --- Verb: is ---
        (1, 3, 1),  # The conversation is playful. (A valid quality for a conversation)
        (3, 3, 1),  # Our Tuna is playful. (A valid quality for a cat)
    ]
    
    print("The true and meaningful sentences are:")
    print("---------------------------------------")

    counts_per_verb = [0] * len(verbs)
    total_count = 0

    for s_idx, v_idx, o_idx in sorted(true_combinations, key=lambda x: x[1]):
        sentence = f"{subjects[s_idx]} {verbs[v_idx]} {objects[o_idx]}."
        print(sentence)
        counts_per_verb[v_idx] += 1
        total_count += 1
    
    print("\n---------------------------------------")
    print("Calculation of the total count:")
    
    equation_parts = []
    for verb, count in zip(verbs, counts_per_verb):
        if count > 0:
            equation_parts.append(str(count))

    equation_str = " + ".join(equation_parts)
    print(f"Number of true sentences = {equation_str} = {total_count}")

solve_sentence_puzzle()
<<<11>>>
def solve_sentence_puzzle():
    """
    This function solves the word puzzle by programmatically identifying all
    true and meaningful sentences that can be formed by swapping subjects and
    objects from the initial set.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    # The verbs are tied to the original sentence structures.
    # To simplify, we'll align verbs with subjects/objects,
    # treating the sentence frame [Subject] [Verb] [Object] as the unit of permutation.
    # Verb index will map to the verb from the original sentence structure.
    # 0: baptized, 1: turned, 2: contains, 3: is
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # This set represents our knowledge base of "true and meaningful" combinations.
    # Each tuple is a (subject_index, verb_index, object_index) combination.
    # The verb_index corresponds to the verb from the original sentence structure.
    # e.g., (0, 3, 3) means "Satyendranath Bose" (subj 0) + "is" (verb from sent 4) + "a Bengalese" (obj 3)
    true_combinations = {
        # Original four true sentences
        (0, 0, 0),  # Satyendranath Bose baptized Bosons.
        (1, 1, 1),  # The conversation turned playful.
        (2, 2, 2),  # Sushi contains fish.
        (3, 3, 3),  # Our Tuna is a Bengalese. (Context: Tuna is a Bengal cat)

        # Newly formed true and meaningful sentences
        (0, 3, 3),  # Satyendranath Bose is a Bengalese. (Factually true)
        (1, 3, 1),  # The conversation is playful. (Synonymous with the original)
        (2, 3, 2),  # Sushi is fish. (Generally true statement)
        (3, 1, 1),  # Our Tuna turned playful. (Plausible for a cat)
        (3, 3, 1),  # Our Tuna is playful. (A cat can be playful)
    }

    found_sentences = []
    # Iterate through all possible sentence structures (defined by the verb)
    for verb_idx, verb in enumerate(verbs):
        # Iterate through all subjects
        for subject_idx, subject in enumerate(subjects):
            # Iterate through all objects
            for object_idx, obj in enumerate(objects):
                if (subject_idx, verb_idx, object_idx) in true_combinations:
                    sentence = f"{subject} {verb} {obj}."
                    found_sentences.append(sentence)
    
    # Sort for consistent output
    found_sentences.sort()

    print("The true and meaningful sentences are:")
    for sentence in found_sentences:
        print(f"- {sentence}")

    print("\nCalculating the total count:")
    
    count_list = [1] * len(found_sentences)
    equation = " + ".join(map(str, count_list))
    
    print(f"{equation} = {len(found_sentences)}")

solve_sentence_puzzle()
<<<9>>>
def solve_sentence_puzzle():
    """
    This function generates all possible sentences by swapping subjects and objects
    from a given set of four sentences, and then identifies which of them are
    true and meaningful based on real-world knowledge.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # Based on real-world knowledge and context, this set contains all the
    # combinations that are considered "true and meaningful".
    true_and_meaningful_sentences = {
        # Group 1: Sentences using the verb 'baptized'
        "Satyendranath Bose baptized Bosons",
        # Group 2: Sentences using the verb 'turned'
        "The conversation turned playful",
        "Our Tuna turned playful",
        # Group 3: Sentences using the verb 'contains'
        "Sushi contains fish",
        "Our Tuna contains fish",
        # Group 4: Sentences using the verb 'is'
        "Satyendranath Bose is a Bengalese",
        "The conversation is playful",
        "Our Tuna is playful",
        "Our Tuna is a Bengalese",
    }

    counts_per_verb = []
    
    print("The true and meaningful sentences are:")
    # Iterate through each verb to represent the four original sentence structures
    for verb in verbs:
        count_for_this_verb = 0
        # For each structure, generate all 16 possible sentences
        for subject in subjects:
            for obj in objects:
                sentence = f"{subject} {verb} {obj}"
                if sentence in true_and_meaningful_sentences:
                    print(f"- {sentence}")
                    count_for_this_verb += 1
        counts_per_verb.append(count_for_this_verb)

    # Calculate the total and format the final equation string
    total_count = sum(counts_per_verb)
    equation_str = " + ".join(map(str, counts_per_verb))

    print("\nFinal Calculation:")
    print(f"The number of true sentences for each verb group are {', '.join(map(str, counts_per_verb))}.")
    print(f"The total is calculated as: {equation_str} = {total_count}")


solve_sentence_puzzle()
<<<9>>>
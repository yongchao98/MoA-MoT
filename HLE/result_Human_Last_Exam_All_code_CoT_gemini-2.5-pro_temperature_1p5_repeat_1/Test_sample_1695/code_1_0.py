def solve_sentence_puzzle():
    """
    Solves the sentence combination puzzle by building all possible sentences
    and checking them against a predefined set of true and meaningful ones.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verb_phrases = ["baptized", "turned", "contains", "is a"]
    objects = ["Bosons", "playful", "fish", "Bengalese"]

    # Based on real-world knowledge and the context provided, these are the
    # sentences determined to be "true and meaningful".
    true_sentences_set = {
        "Satyendranath Bose baptized Bosons.",
        "The conversation turned playful.",
        "Sushi contains fish.",
        "Our Tuna is a Bengalese.",
        "Satyendranath Bose is a Bengalese.",
        "Our Tuna turned playful.",
        "Our Tuna contains fish.",
    }

    count = 0
    equation_parts = []

    print("The true and meaningful sentences are:")
    print("------------------------------------")

    # Iterate through all combinations to find the true ones
    for s in subjects:
        for v in verb_phrases:
            for o in objects:
                # Construct the sentence
                sentence = f"{s} {v} {o}."
                
                # Check if the constructed sentence is in our set of true sentences
                if sentence in true_sentences_set:
                    count += 1
                    equation_parts.append("1")
                    print(sentence)

    # Print the final calculation and result
    print("\n------------------------------------")
    print("Calculating the total count:")
    
    if count > 0:
        equation_string = " + ".join(equation_parts)
        print(f"{equation_string} = {count}")
    else:
        print("No true and meaningful sentences were found.")
        print(f"Total count = {count}")

solve_sentence_puzzle()

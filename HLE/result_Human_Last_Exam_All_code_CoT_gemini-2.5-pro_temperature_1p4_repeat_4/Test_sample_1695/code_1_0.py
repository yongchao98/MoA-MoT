def solve_sentence_puzzle():
    """
    Solves the sentence logic puzzle by identifying all true and meaningful
    combinations from a given set of subjects, verbs, and objects.
    """

    # 1. Deconstruct the sentences into their components
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # 2. Define the set of sentences that are true and meaningful
    # This is based on logical deduction and real-world knowledge.
    true_and_meaningful_sentences = {
        "Satyendranath Bose baptized Bosons.",
        "The conversation turned playful.",
        "Sushi contains fish.",
        "Our Tuna is a Bengalese.",
        "The conversation is playful.",
        "Our Tuna is playful.",
        "Our Tuna turned playful."
    }

    found_sentences = []
    
    print("The true and meaningful sentences are:")
    # 3. Iterate through all possible combinations and check against our set
    for s in subjects:
        for v in verbs:
            for o in objects:
                # Form the sentence
                sentence = f"{s} {v} {o}."
                
                # Check if it's in our pre-determined list of true sentences
                if sentence in true_and_meaningful_sentences:
                    found_sentences.append(sentence)
                    
    # Print the sentences found and the final equation
    for sentence in sorted(found_sentences): # Sorting for consistent output
        print(sentence)
        
    count = len(found_sentences)
    equation_parts = ["1"] * count
    
    print("\nTotal count calculation:")
    print(f"{' + '.join(equation_parts)} = {count}")


solve_sentence_puzzle()
<<<7>>>
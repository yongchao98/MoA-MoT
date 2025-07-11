def solve_sentence_puzzle():
    """
    Solves the sentence swapping puzzle by generating all possible sentences
    and checking them against a pre-determined set of true and meaningful ones.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]
    verb_structures = ["baptized", "turned", "contains", "is"]

    # Based on world knowledge, these are the true and meaningful sentences that can be formed.
    true_sentences = {
        "Satyendranath Bose baptized Bosons.",
        "The conversation turned playful.",
        "Our Tuna turned playful.",
        "Sushi contains fish.",
        "Satyendranath Bose is a Bengalese.",
        "The conversation is playful.",
        "Our Tuna is playful.",
        "Our Tuna is a Bengalese.",
    }

    print("The true and meaningful sentences are:\n")

    counts_per_structure = []
    
    for i, verb in enumerate(verb_structures):
        structure_count = 0
        
        # This block is just for printing a clean header for each group
        if i > 0:
            print()

        for s in subjects:
            for o in objects:
                # Generate the sentence with basic grammar rules
                sentence = ""
                if verb == "is":
                    # Adjectives don't take an article 'a'
                    if o == "playful":
                        sentence = f"{s} is {o}."
                    else:
                        sentence = f"{s} is {o}."
                else:
                    sentence = f"{s} {verb} {o}."

                if sentence in true_sentences:
                    print(sentence)
                    structure_count += 1
        
        counts_per_structure.append(structure_count)

    total_count = sum(counts_per_structure)

    # Format the final equation string
    equation_str = " + ".join(map(str, counts_per_structure))

    print("\n------------------------------------")
    print(f"Calculation: {equation_str} = {total_count}")
    print(f"Total number of true and meaningful sentences: {total_count}")

solve_sentence_puzzle()
<<<8>>>
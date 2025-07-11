def solve_sentence_puzzle():
    """
    This function generates all possible sentences by swapping subjects and objects,
    evaluates which ones are true and meaningful, and prints the result.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # Based on logical analysis, these are the only true and meaningful combinations.
    # The originals are included.
    known_true_sentences = {
        "Satyendranath Bose baptized Bosons.",
        "The conversation turned playful.",
        "Sushi contains fish.",
        "Our Tuna is a Bengalese.",
        "Our Tuna is playful." # This is the only new combination that is true.
    }

    found_true_sentences = []
    count_components = []

    # Generate all 4x4 = 16 sentences
    for i in range(len(subjects)):
        for j in range(len(objects)):
            # The verb is associated with its original subject
            sentence = f"{subjects[i]} {verbs[i]} {objects[j]}."
            
            if sentence in known_true_sentences:
                found_true_sentences.append(sentence)
                count_components.append("1")

    print("The true and meaningful sentences are:")
    for i, sentence in enumerate(found_true_sentences, 1):
        print(f"{i}. {sentence}")

    equation = " + ".join(count_components)
    total = len(found_true_sentences)
    
    print(f"\nTotal count: {equation} = {total}")

solve_sentence_puzzle()
<<<5>>>
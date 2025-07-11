def solve_sentence_puzzle():
    """
    Solves the sentence swapping puzzle by evaluating all possible combinations.
    """
    subjects = [
        "Satyendranath Bose",
        "The conversation",
        "Sushi",
        "Our Tuna"
    ]

    predicates = [
        "baptized Bosons.",
        "turned playful.",
        "contains fish.",
        "is a Bengalese."
    ]

    # A matrix to represent the truth value of each combination (subject, predicate)
    # Rows correspond to subjects, columns to predicates.
    # True means the sentence is true and meaningful.
    evaluations = [
        # Predicates: baptized Bosons, turned playful, contains fish, is a Bengalese
        [True, False, False, True],   # Subject: Satyendranath Bose
        [False, True, False, False],  # Subject: The conversation
        [False, False, True, False],  # Subject: Sushi
        [False, True, True, True]     # Subject: Our Tuna
    ]

    true_sentences = []
    for i in range(len(subjects)):
        for j in range(len(predicates)):
            if evaluations[i][j]:
                sentence = f"{subjects[i]} {predicates[j]}"
                true_sentences.append(sentence)

    print("The true and meaningful sentences are:")
    for sentence in true_sentences:
        print(f"- {sentence}")

    # Build and print the final equation as requested
    count = len(true_sentences)
    equation_parts = ["1"] * count
    equation = " + ".join(equation_parts)
    print(f"\nThe total count is {equation} = {count}")


solve_sentence_puzzle()
print("\n<<<7>>>")
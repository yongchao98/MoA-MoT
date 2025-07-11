def solve_sentence_puzzle():
    """
    This function forms new sentences by combining subjects, verbs, and objects
    from a given set of four sentences and counts how many are true and meaningful.
    """
    # Step 1: Define the components of the original sentences.
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # Step 2: Identify the combinations (subject_index, verb_index, object_index)
    # that form true and meaningful sentences based on real-world and contextual knowledge.
    # Indices are 0-based.
    true_sentence_indices = [
        (0, 0, 0),  # Satyendranath Bose baptized Bosons. (Original)
        (1, 1, 1),  # The conversation turned playful. (Original)
        (2, 2, 2),  # Sushi contains fish. (Original)
        (3, 3, 3),  # Our Tuna is a Bengalese. (Original)
        (3, 1, 1),  # Our Tuna turned playful. (New, but true for a cat)
        (3, 3, 1)   # Our Tuna is playful. (New, but true for a cat)
    ]

    print("Found the following true and meaningful sentences:")
    
    count = 0
    # Step 3: Iterate through the identified true sentences and print them.
    for s_idx, v_idx, o_idx in true_sentence_indices:
        count += 1
        sentence = f"{subjects[s_idx]} {verbs[v_idx]} {objects[o_idx]}."
        # The line below shows the "equation" for each found sentence
        print(f"Sentence {count}: {sentence}")

    print("\n-------------------------------------------------")
    print(f"The total number of true and meaningful sentences is the sum of each sentence found:")
    equation = " + ".join(["1"] * count)
    print(f"{equation} = {count}")
    print("-------------------------------------------------")


solve_sentence_puzzle()
<<<6>>>
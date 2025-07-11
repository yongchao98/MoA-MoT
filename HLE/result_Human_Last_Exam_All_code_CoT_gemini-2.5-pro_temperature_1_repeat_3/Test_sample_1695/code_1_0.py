def solve_sentence_puzzle():
    """
    This function identifies and counts the number of true and meaningful sentences
    that can be formed by combining subjects, verbs, and objects from a given set.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # Triples of (subject_index, verb_index, object_index) for sentences
    # that are determined to be true and meaningful.
    true_sentence_indices = [
        # 1. The original sentence: "Satyendranath Bose baptized Bosons."
        # This is factually true in a figurative sense (he is the namesake of bosons).
        (0, 0, 0),

        # 2. The original sentence: "The conversation turned playful."
        # True by premise.
        (1, 1, 1),

        # 3. The original sentence: "Sushi contains fish."
        # True by general knowledge.
        (2, 2, 2),

        # 4. The original sentence: "Our Tuna is a Bengalese."
        # True by premise, establishing "Tuna" as a cat of the Bengalese breed.
        (3, 3, 3),

        # 5. New sentence: "The conversation is playful."
        # "is" and "turned" are similar linking verbs; this is a true rephrasing.
        (1, 3, 1),

        # 6. New sentence: "Our Tuna turned playful."
        # We know Tuna is a cat, and it is a general truth that cats can be playful.
        (3, 1, 1),
        
        # 7. New sentence: "Our Tuna is playful."
        # Similar to the above, this is a statement about a cat's trait.
        (3, 3, 1),

        # 8. New sentence: "Our Tuna contains fish."
        # The cat's name "Tuna" strongly implies its diet includes fish.
        (3, 2, 2),
    ]

    print("The true and meaningful sentences are:")
    print("-" * 35)

    equation_parts = []
    for idx_tuple in true_sentence_indices:
        s_idx, v_idx, o_idx = idx_tuple
        sentence = f"{subjects[s_idx]} {verbs[v_idx]} {objects[o_idx]}."
        print(sentence)
        equation_parts.append("1")
    
    print("-" * 35)
    total_count = len(equation_parts)
    equation_str = " + ".join(equation_parts)
    print(f"Final calculation: {equation_str} = {total_count}")


solve_sentence_puzzle()

<<<8>>>
def solve_sentence_puzzle():
    """
    Solves the sentence puzzle by identifying and counting true and meaningful
    sentences formed by swapping subjects and objects.
    """
    # Step 1: Define the components of the sentences
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]
    
    # Sentence structures (verbs/predicates)
    structures = [
        "baptized",
        "turned",
        "contains",
        "is"
    ]

    # Step 2: Manually determine which combinations are "true and meaningful"
    # This list of tuples represents the valid (subject_index, structure_index, object_index)
    # based on world knowledge and logic.
    true_combinations = [
        # Structure 1: '... baptized ...'
        (0, 0, 0),  # Satyendranath Bose baptized Bosons. (True, namesake)
        
        # Structure 2: '... turned ...'
        (1, 1, 1),  # The conversation turned playful. (True, idiomatic)
        (3, 1, 1),  # Our Tuna turned playful. (True, plausible for a cat)
        
        # Structure 3: '... contains ...'
        (2, 2, 2),  # Sushi contains fish. (True, common ingredient)
        (3, 2, 2),  # Our Tuna contains fish. (True, plausible for a cat's diet)

        # Structure 4: '... is ...'
        (1, 3, 1),  # The conversation is playful. (True, descriptive)
        (2, 3, 2),  # Sushi is fish. (True, common description)
        (3, 3, 1),  # Our Tuna is playful. (True, plausible for a cat)
        (3, 3, 3)   # Our Tuna is a Bengalese. (True, statement of breed)
    ]

    print("The true and meaningful sentences are:")
    print("-" * 35)

    counts_per_structure = {s: 0 for s in structures}
    
    # Step 3: Print the valid sentences and count them
    for s_idx, v_idx, o_idx in true_combinations:
        subject = subjects[s_idx]
        structure = structures[v_idx]
        obj = objects[o_idx]
        
        # Print the reconstructed sentence
        print(f"{subject} {structure} {obj}.")
        
        # Increment the count for the respective structure
        counts_per_structure[structure] += 1

    print("-" * 35)
    
    # Step 4: Display the final calculation
    # "Remember in the final code you still need to output each number in the final equation!"
    
    counts = list(counts_per_structure.values())
    total = sum(counts)
    
    # Create the equation string
    equation = " + ".join(map(str, counts))
    
    print(f"Number of true sentences per structure:")
    for i, structure in enumerate(structures):
        print(f"  - Sentences with '{structure}': {counts[i]}")

    print(f"\nTotal calculation: {equation} = {total}")

solve_sentence_puzzle()
<<<9>>>
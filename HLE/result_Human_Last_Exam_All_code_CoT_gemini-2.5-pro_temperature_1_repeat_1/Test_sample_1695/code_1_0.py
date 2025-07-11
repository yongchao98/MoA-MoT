def solve_sentence_puzzle():
    """
    Generates all possible sentences by swapping subjects and objects from the
    four base sentences and counts how many are true and meaningful.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]
    
    # The predicates define the four sentence structures.
    predicates = ["baptized", "turned", "contains", "is a"]
    
    true_sentences = []
    # A dictionary to store the count of true sentences for each predicate.
    predicate_counts = {p: 0 for p in predicates}

    # Systematically check all 4x4x4 = 64 possible sentences.
    for predicate in predicates:
        for subject in subjects:
            for obj in objects:
                is_true = False
                
                # Rule 1: 'baptized' - only Bose naming Bosons is true.
                if predicate == "baptized":
                    if subject == "Satyendranath Bose" and obj == "Bosons":
                        is_true = True
                
                # Rule 2: 'turned' - a person, conversation, or cat can become playful.
                elif predicate == "turned":
                    if subject in ["Satyendranath Bose", "The conversation", "Our Tuna"] and obj == "playful":
                        is_true = True
                
                # Rule 3: 'contains' - only sushi contains fish.
                elif predicate == "contains":
                    if subject == "Sushi" and obj == "fish":
                        is_true = True
                
                # Rule 4: 'is a' - "Our Tuna" is established as a Bengalese cat.
                elif predicate == "is a":
                    if subject == "Our Tuna" and obj == "a Bengalese":
                        is_true = True
                
                if is_true:
                    sentence = f"{subject} {predicate} {obj}."
                    true_sentences.append(sentence)
                    predicate_counts[predicate] += 1

    print("The true and meaningful sentences are:")
    for sentence in true_sentences:
        print(sentence)
    print("\n---")
    
    # Retrieve the counts in the order of the predicates list for the final equation.
    counts = [predicate_counts[p] for p in predicates]
    total = sum(counts)
    
    # Create the equation string as requested.
    equation_parts = [str(c) for c in counts]
    equation = " + ".join(equation_parts)

    print(f"The calculation is based on the number of true sentences found for each of the 4 structures:")
    print(f"Structure 1 ('baptized'): {counts[0]}")
    print(f"Structure 2 ('turned'):   {counts[1]}")
    print(f"Structure 3 ('contains'): {counts[2]}")
    print(f"Structure 4 ('is a'):     {counts[3]}")
    print("\nFinal Equation:")
    print(f"{equation} = {total}")

    print(f"\n<<<{total}>>>")

solve_sentence_puzzle()
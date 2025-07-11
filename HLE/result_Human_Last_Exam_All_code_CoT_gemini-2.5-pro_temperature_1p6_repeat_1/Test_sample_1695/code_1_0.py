def solve_sentence_puzzle():
    """
    Generates all possible sentences by swapping subjects and objects from four
    base sentences, and counts how many are true and meaningful.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is a"]
    objects = ["Bosons", "playful", "fish", "Bengalese"]

    true_sentences = []
    # Using a dictionary to store counts for each verb category
    counts = {verb: 0 for verb in verbs}

    # Iterate through all combinations
    for subject in subjects:
        for verb in verbs:
            for obj in objects:
                is_true = False
                # Apply logical rules based on our model of the world
                if verb == "baptized":
                    # Only a person can baptize (name) something, and we know S.Bose "baptized" Bosons.
                    if subject == "Satyendranath Bose" and obj == "Bosons":
                        is_true = True
                
                elif verb == "turned":
                    # Something 'turns' into a state/quality.
                    if obj == "playful":
                        # A conversation or a pet can turn playful.
                        if subject in ["The conversation", "Our Tuna"]:
                            is_true = True

                elif verb == "contains":
                    # Sushi contains fish. A cat (Tuna) can contain fish (if it ate it).
                    if obj == "fish":
                        if subject in ["Sushi", "Our Tuna"]:
                            is_true = True
                            
                elif verb == "is a":
                    # Checking identity or category.
                    # S. Bose was a person of Bengali ethnicity.
                    if subject == "Satyendranath Bose" and obj == "Bengalese":
                        is_true = True
                    # Our pet cat Tuna is of the Bengalese breed.
                    if subject == "Our Tuna" and obj == "Bengalese":
                        is_true = True

                if is_true:
                    sentence = f"{subject} {verb} {obj}."
                    true_sentences.append(sentence)
                    counts[verb] += 1

    print("The true and meaningful sentences are:")
    for sentence in sorted(true_sentences):
        print(f"- {sentence}")
    
    print("\n" + "="*40)
    print("Final Count Calculation:")
    print("="*40)

    total_count = sum(counts.values())
    calculation_str = " + ".join(str(c) for c in counts.values())

    print(f"Number of valid sentences using 'baptized': {counts['baptized']}")
    print(f"Number of valid sentences using 'turned':   {counts['turned']}")
    print(f"Number of valid sentences using 'contains': {counts['contains']}")
    print(f"Number of valid sentences using 'is a':     {counts['is a']}")
    
    print("\nThe total is the sum of these counts:")
    # The final equation as requested
    print(f"{counts['baptized']} + {counts['turned']} + {counts['contains']} + {counts['is a']} = {total_count}")


solve_sentence_puzzle()

<<<7>>>
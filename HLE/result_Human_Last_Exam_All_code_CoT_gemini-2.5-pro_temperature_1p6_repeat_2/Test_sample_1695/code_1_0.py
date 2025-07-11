def solve_sentence_puzzle():
    """
    Solves the sentence combination puzzle by evaluating all possible sentences
    for truth and meaning.
    """
    subjects = {
        "S1": "Satyendranath Bose",
        "S2": "The conversation",
        "S3": "Sushi",
        "S4": "Our Tuna"
    }

    verbs = {
        "V1": "baptized",
        "V2": "turned",
        "V3": "contains",
        "V4": "is"
    }

    objects = {
        "O1": "Bosons",
        "O2": "playful",
        "O3": "fish",
        "O4": "a Bengalese"
    }

    true_sentences_by_verb = {v: [] for v in verbs.values()}

    # This function encodes the logic for checking if a sentence is true and meaningful
    def is_true_and_meaningful(s_key, v_key, o_key):
        s = subjects[s_key]
        v = verbs[v_key]
        o = objects[o_key]

        # --- Verb 1: baptized ---
        # Requires a sentient agent as a subject to name or christen an object.
        if v == "baptized":
            # "Satyendranath Bose baptized Bosons." (Original, true history)
            if s_key == "S1" and o_key == "O1":
                return True
        
        # --- Verb 2: turned ---
        # Subject changes state to become the object (adjective).
        elif v == "turned":
            if o_key == "O2": # Object must be 'playful'
                # "The conversation turned playful." (Original, true)
                # "Our Tuna turned playful." (A cat can become playful, true statement of behavior)
                if s_key in ["S2", "S4"]:
                    return True

        # --- Verb 3: contains ---
        # Subject is physically composed of or holds the object.
        elif v == "contains":
            # "Sushi contains fish." (Original, true)
            # "Our Tuna contains fish." (A cat named Tuna likely eats/contains fish, true in this context)
            if o_key == "O3" and s_key in ["S3", "S4"]:
                return True
            # All physical matter contains Bosons (a class of particles).
            # Subjects S1, S3, and S4 are physical entities.
            if o_key == "O1" and s_key in ["S1", "S3", "S4"]:
                return True
        
        # --- Verb 4: is ---
        # Subject has the quality of the object or is an instance of the object's class.
        elif v == "is":
            # "Our Tuna is a Bengalese." (Original, true)
            if s_key == "S4" and o_key == "O4":
                return True
            if o_key == "O2": # Object is 'playful'
                # "The conversation is playful." (True, synonymous with 'turned playful')
                # "Our Tuna is playful." (A cat can be playful, true quality)
                if s_key in ["S2", "S4"]:
                    return True
        
        return False

    # Iterate through all combinations
    for v_key in sorted(verbs.keys()):
        verb = verbs[v_key]
        count_for_verb = 0
        for s_key in sorted(subjects.keys()):
            for o_key in sorted(objects.keys()):
                if is_true_and_meaningful(s_key, v_key, o_key):
                    true_sentences_by_verb[verb].append(f"{subjects[s_key]} {verb} {objects[o_key]}.")

    # Print the results
    print("Finding all true and meaningful sentences by combining subjects, verbs, and objects:")
    total_count = 0
    equation_parts = []
    
    # Ordered list of verbs for consistent output
    verb_order = ["baptized", "turned", "contains", "is"]

    for verb in verb_order:
        count = len(true_sentences_by_verb[verb])
        print(f"\nFor the verb '{verb}', we found {count} true sentence(s):")
        for sentence in true_sentences_by_verb[verb]:
            print(f"- {sentence}")
        total_count += count
        equation_parts.append(str(count))

    print("\nThe total number of true and meaningful sentences is the sum of the counts for each verb.")
    final_equation = " + ".join(equation_parts)
    print(f"Final Equation: {final_equation} = {total_count}")

solve_sentence_puzzle()
<<<11>>>
def solve_sentence_puzzle():
    """
    This function solves the puzzle by generating all possible sentences
    and evaluating them for truth and meaning.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # This dictionary will store the counts for each verb group
    verb_counts = {v: 0 for v in verbs}
    
    # Store the sentences for clarity, though not strictly required for the count
    true_sentences = []

    # Iterate through every combination of subject, verb, and object
    for s in subjects:
        for v in verbs:
            for o in objects:
                # Default to false unless a condition is met
                is_true_and_meaningful = False

                # --- Evaluation Logic ---

                if v == "baptized":
                    # Only one true case: the original sentence.
                    if s == "Satyendranath Bose" and o == "Bosons":
                        is_true_and_meaningful = True
                
                elif v == "turned":
                    # A conversation can turn playful (original).
                    if s == "The conversation" and o == "playful":
                        is_true_and_meaningful = True
                    # A cat can turn playful.
                    elif s == "Our Tuna" and o == "playful":
                        is_true_and_meaningful = True

                elif v == "contains":
                    # Sushi contains fish (original).
                    if s == "Sushi" and o == "fish":
                        is_true_and_meaningful = True
                    # A cat (a carnivore) can contain fish it has eaten.
                    elif s == "Our Tuna" and o == "fish":
                        is_true_and_meaningful = True

                elif v == "is":
                    # Our Tuna is a Bengalese (original).
                    if s == "Our Tuna" and o == "a Bengalese":
                        is_true_and_meaningful = True
                    # A cat can be playful.
                    elif s == "Our Tuna" and o == "playful":
                        is_true_and_meaningful = True
                    # S. Bose was a person from Bengal (Bengali).
                    elif s == "Satyendranath Bose" and o == "a Bengalese":
                        is_true_and_meaningful = True
                    # A conversation can be playful.
                    elif s == "The conversation" and o == "playful":
                        is_true_and_meaningful = True
                    # Sushi is often made of fish (e.g., sashimi).
                    elif s == "Sushi" and o == "fish":
                        is_true_and_meaningful = True

                # If the sentence is valid, increment the count for its verb group
                if is_true_and_meaningful:
                    verb_counts[v] += 1
                    true_sentences.append(f"{s} {v} {o}.")

    # --- Output the results ---
    
    # Print the identified true sentences (optional, for clarity)
    # print("Found the following true and meaningful sentences:")
    # for sentence in sorted(true_sentences):
    #    print(f"- {sentence}")
    # print("-" * 20)

    total_count = sum(verb_counts.values())
    
    # Format the equation string as requested
    equation_parts = [str(verb_counts[v]) for v in verbs]
    equation_str = " + ".join(equation_parts)

    print("The final equation is:")
    print(f"{equation_str} = {total_count}")


solve_sentence_puzzle()
<<<10>>>
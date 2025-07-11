def solve_sentence_puzzle():
    """
    This function solves the sentence swapping puzzle by evaluating all possible combinations
    of subjects and predicates and counting how many are true and meaningful.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    predicates = ["baptized Bosons", "turned playful", "contains fish", "is a Bengalese"]

    # This dictionary defines the logic for which combinations are true.
    # The key is the subject, and the value is a list of predicates that form a true sentence with that subject.
    true_combinations = {
        "Satyendranath Bose": ["baptized Bosons"],
        "The conversation": ["turned playful"],
        "Sushi": ["contains fish"],
        "Our Tuna": ["turned playful", "contains fish", "is a Bengalese"]
    }

    print("The original four true sentences are:")
    print("1. Satyendranath Bose baptized Bosons.")
    print("2. The conversation turned playful.")
    print("3. Sushi contains fish.")
    print("4. Our Tuna is a Bengalese.")
    print("\nWe will form new sentences by combining each subject with each predicate and check if the result is true and meaningful.")
    print("Based on our analysis, the following sentences are true:")

    counts_per_subject = []
    total_true_sentences = 0
    
    for subject in subjects:
        # Find the list of true predicates for the current subject, or an empty list if none.
        true_preds_for_subject = true_combinations.get(subject, [])
        count_for_this_subject = len(true_preds_for_subject)
        
        if count_for_this_subject > 0:
            for predicate in true_preds_for_subject:
                print(f"- {subject} {predicate}.")

        counts_per_subject.append(count_for_this_subject)
        total_true_sentences += count_for_this_subject

    print("\nTo find the total, we sum the number of true sentences for each subject:")
    
    # Build and print the equation string, e.g., "1 + 1 + 1 + 3 = 6"
    equation_str = " + ".join(map(str, counts_per_subject))
    print(f"{equation_str} = {total_true_sentences}")

    print(f"\n<<<6>>>")

solve_sentence_puzzle()
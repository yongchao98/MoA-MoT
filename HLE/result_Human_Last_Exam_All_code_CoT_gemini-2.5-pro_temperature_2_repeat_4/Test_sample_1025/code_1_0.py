def rerank_documents():
    """
    This function analyzes the relevance of two documents to two questions and determines the best ranking.
    """

    # --- Analysis ---
    # Question 1: "What animal's coat color changed as a baby?"
    # The key distinction between the documents is the timing of the coat color change for young bulls.
    # Doc 1: "at a few months old" - This aligns well with the term "baby".
    # Doc 2: "as they reached puberty" - This refers to a later, adolescent stage, not a baby stage.
    # Conclusion for Q1: Document 1 is more relevant.
    rank_for_q1 = 1

    # Question 2: "What animal's coat color did not change as a baby?"
    # Both documents state that "Cows retained the reddish-brown colour."
    # However, Doc 1 provides a better contextual timeline ("a few months old") that is relevant to the "baby" stage mentioned in the question.
    # Doc 2's timeline ("puberty") is less relevant to the question's focus.
    # Conclusion for Q2: Document 1 provides a better, more relevant context.
    rank_for_q2 = 1

    # --- Output ---
    # The final answer format combines the ranks for Q1 and Q2.
    final_answer = f"{rank_for_q1}-{rank_for_q2}"

    print("This script determines the best document for each question.")
    print("Based on the analysis, here is the resulting 'equation':")
    print(f"Rank for Question 1: {rank_for_q1}")
    print(f"Rank for Question 2: {rank_for_q2}")
    print(f"Final Combined Answer: {final_answer}")


rerank_documents()
<<<A>>>
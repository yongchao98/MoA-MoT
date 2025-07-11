def solve_reranking():
    """
    This function analyzes the relevance of two documents for two questions and determines the final ranking.
    """

    # Analysis for Question 1: "What animal's coat color changed as a baby?"
    q1_doc1_evidence = "Young bulls changed their coat colour at a few months old..."
    q1_doc2_evidence = "As they reached puberty, young bulls changed their coat colour..."
    # "at a few months old" is a more direct and accurate answer for "as a baby" than "As they reached puberty".
    # Therefore, Document 1 is more relevant for Question 1.
    q1_rank = 1

    # Analysis for Question 2: "What animal's coat color did not change as a baby?"
    q2_doc1_evidence = "Cows retained the reddish-brown colour."
    q2_doc2_evidence = "Cows retained the reddish-brown colour."
    # Both documents provide the exact same sentence, which perfectly answers the question. They are equally relevant.
    # By standard tie-breaking rules in information retrieval, when relevance is equal, the first-listed document is chosen.
    # Therefore, we choose Document 1 for Question 2.
    q2_rank = 1

    # The final combined ranking is "1-1". This corresponds to answer choice A.
    
    print("The ranking for Question 1 is: " + str(q1_rank))
    print("The ranking for Question 2 is: " + str(q2_rank))
    print("The final combined ranking is an equation showing the rank for each question.")
    print(f"Final Equation: {q1_rank}-{q2_rank}")

solve_reranking()
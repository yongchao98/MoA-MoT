def solve_reranking():
    """
    This function analyzes two documents to determine the best fit for two questions,
    simulating an information retrieval reranking task.
    """

    # Define the documents
    doc1_text = "Young bulls changed their coat colour at a few months old... Cows retained the reddish-brown colour."
    doc2_text = "As they reached puberty, young bulls changed their coat colour... Cows retained the reddish-brown colour."

    # --- Reranking for Question 1: "What animal's coat color changed as a baby?" ---
    # We score based on the relevance of the time frame to "as a baby".
    # "a few months old" is a much better match than "puberty".
    score1_q1 = 2 if "a few months old" in doc1_text else 0
    score2_q1 = 1 if "puberty" in doc2_text else 0

    if score1_q1 > score2_q1:
        rank_q1 = 1
    elif score2_q1 > score1_q1:
        rank_q1 = 2
    else:
        # Default tie-breaker (prefer lower doc number)
        rank_q1 = 1

    # --- Reranking for Question 2: "What animal's coat color did not change as a baby?" ---
    # We score based on the presence of the keyword "retained".
    score1_q2 = 1 if "retained" in doc1_text else 0
    score2_q2 = 1 if "retained" in doc2_text else 0
    
    if score1_q2 > score2_q2:
        rank_q2 = 1
    elif score2_q2 > score1_q2:
        rank_q2 = 2
    else:
        # Tie-breaker: Since scores are equal, we favor the document
        # that ranked higher for the first question.
        rank_q2 = rank_q1

    # --- Output the final answer ---
    # The final answer format is "rank_for_Q1-rank_for_Q2".
    # In our case, rank_for_Q1 is 1 and rank_for_Q2 is 1.
    print("The final ranking is an equation of 'Rank for Question 1 - Rank for Question 2'")
    print(f"The number for the first part of the equation is: {rank_q1}")
    print(f"The number for the second part of the equation is: {rank_q2}")
    print(f"Final Equation: {rank_q1}-{rank_q2}")

solve_reranking()
<<<A>>>
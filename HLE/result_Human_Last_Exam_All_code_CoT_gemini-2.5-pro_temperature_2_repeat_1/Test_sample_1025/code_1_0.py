def solve_reranking():
    """
    Solves the document reranking task by analyzing relevance for each question.
    """
    # Plan:
    # 1. Analyze Question 1: Determine which document is more relevant.
    #    - Question 1: "What animal's coat color changed as a baby?"
    #    - Doc 1 states the change happened "at a few months old".
    #    - Doc 2 states the change happened "As they reached puberty".
    #    - "a few months old" is a better match for "baby" than "puberty".
    #    - Therefore, Doc 1 is top-ranked for Q1.
    q1_rank = 1

    # 2. Analyze Question 2: Determine which document is more relevant.
    #    - Question 2: "What animal's coat color did not change as a baby?"
    #    - Both Doc 1 and Doc 2 contain the identical sentence: "Cows retained the reddish-brown colour."
    #    - This makes them equally relevant. As a tie-breaker, we can prefer the document
    #      that was ranked higher for the other question.
    #    - Therefore, Doc 1 is chosen as the top-ranked document for Q2.
    q2_rank = 1

    # 3. Combine and print the results.
    print("Ranking for Question 1: The phrase 'at a few months old' in Document 1 is more relevant to 'as a baby' than 'puberty' in Document 2.")
    print(f"So the choice for Question 1 is: {q1_rank}")
    print("\nRanking for Question 2: Both documents contain the identical answer. Using the rank from Question 1 as a tie-breaker.")
    print(f"So the choice for Question 2 is: {q2_rank}")

    # The final equation is the combination of the two ranks.
    print(f"\nFinal combined answer: {q1_rank}-{q2_rank}")

solve_reranking()
<<<A>>>
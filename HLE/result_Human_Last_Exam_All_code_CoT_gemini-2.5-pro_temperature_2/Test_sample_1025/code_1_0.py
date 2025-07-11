def solve_reranking():
    """
    Analyzes two documents against two questions to determine the best ranking.
    """
    # Analysis for Question 1: What animal's coat color changed as a baby?
    # Doc 1 says: "Young bulls changed their coat colour at a few months old..."
    # Doc 2 says: "As they reached puberty, young bulls changed their coat colour..."
    # Both documents directly answer that young aurochs bulls changed color. They are equally relevant.
    q1_ranking = "equally_ranked" # Let's denote this as 1 for doc 1 and 2 for doc 2 being equal

    # Analysis for Question 2: What animal's coat color did not change as a baby?
    # Doc 1 says: "Cows retained the reddish-brown colour."
    # Doc 2 says: "Cows retained the reddish-brown colour."
    # Both documents contain the identical sentence answering the question. They are equally relevant.
    q2_ranking = "equally_ranked"

    # Conclusion
    # For Question 1, both documents are equally good.
    # For Question 2, both documents are equally good.
    # Therefore, the overall answer is that both documents are equally ranked.
    final_choice = "E"
    
    print(f"Reasoning:")
    print("For Question 1, both documents state that the young bull's coat changed color. Thus, they are equally relevant.")
    print("For Question 2, both documents state that the cow's coat did not change. Thus, they are equally relevant.")
    print(f"Conclusion: Both documents are equally ranked for both questions.")
    print(f"The correct option is {final_choice}.")

solve_reranking()
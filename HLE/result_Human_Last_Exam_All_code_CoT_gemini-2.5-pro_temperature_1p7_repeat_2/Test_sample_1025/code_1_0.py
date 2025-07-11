def solve_reranking_task():
    """
    This function explains the step-by-step reasoning for the document reranking task.
    """

    print("Step 1: Analyze Question 1: 'What animal's coat color changed as a baby?'")
    print("  - Doc 1 mentions the change happened 'at a few months old'.")
    print("  - Doc 2 mentions the change happened when they 'reached puberty'.")
    print("  - The phrase 'at a few months old' is more relevant to 'as a baby' than 'puberty'.")
    print("  - Therefore, for Question 1, Document 1 is ranked higher.")
    q1_rank = 1
    print(f"  - Ranking for Question 1: {q1_rank}\n")

    print("Step 2: Analyze Question 2: 'What animal's coat color did not change as a baby?'")
    print("  - Both Doc 1 and Doc 2 contain the identical, most relevant sentence: 'Cows retained the reddish-brown colour.'")
    print("  - Since the information is identical, the documents are equally relevant.")
    print("  - In case of a tie, a standard tie-breaking rule is to choose the first document.")
    print("  - Therefore, for Question 2, we select Document 1.")
    q2_rank = 1
    print(f"  - Ranking for Question 2: {q2_rank}\n")

    print("Step 3: Combine the rankings.")
    print(f"The combined ranking is '{q1_rank}-{q2_rank}'.")
    print("This corresponds to answer choice A.")

solve_reranking_task()
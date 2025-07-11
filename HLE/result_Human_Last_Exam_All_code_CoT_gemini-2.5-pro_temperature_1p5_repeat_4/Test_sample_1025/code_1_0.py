def solve_reranking():
    """
    Analyzes two documents to find the best match for two questions and provides a final ranking.
    """
    doc1_text = "Young bulls changed their coat colour at a few months old to black..."
    doc2_text = "As they reached puberty, young bulls changed their coat colour to black..."
    doc_shared_text = "Cows retained the reddish-brown colour."

    # --- Step 1: Analyze Question 1 ---
    # Question 1: What animal's coat color changed as a baby?
    # The phrase "at a few months old" (Doc 1) is more relevant to "as a baby"
    # than "as they reached puberty" (Doc 2).
    print("Evaluating Question 1: 'What animal's coat color changed as a baby?'")
    score_q1_doc1 = 1 if "at a few months old" in doc1_text else 0
    score_q1_doc2 = 1 if "at a few months old" in doc2_text else 0 # This will be 0

    if score_q1_doc1 > score_q1_doc2:
        rank_q1 = 1
        print("Finding: Document 1 is more relevant because 'at a few months old' is a better match for 'as a baby'.")
    else:
        rank_q1 = 2 # Failsafe, not expected
        print("Finding: Document 2 is more relevant.")
    print(f"Chosen document for Question 1: {rank_q1}\n")


    # --- Step 2: Analyze Question 2 ---
    # Question 2: What animal's coat color did not change as a baby?
    # The relevant sentence is identical in both documents.
    print("Evaluating Question 2: 'What animal's coat color did not change as a baby?'")
    score_q2_doc1 = 1 if "retained" in doc_shared_text else 0
    score_q2_doc2 = 1 if "retained" in doc_shared_text else 0

    if score_q2_doc1 == score_q2_doc2:
        print("Finding: Both documents contain the same relevant information and are equally ranked.")
        # --- Step 3: Apply Tie-Breaker ---
        # Prefer the document that won the first query.
        rank_q2 = rank_q1
        print(f"Tie-breaking: Since Document {rank_q1} was better for Question 1, we choose it for Question 2 as well.")
    else:
        rank_q2 = 1 if score_q2_doc1 > score_q2_doc2 else 2 # Failsafe
    print(f"Chosen document for Question 2: {rank_q2}\n")


    # --- Step 4: Final Answer ---
    print("The final ranking combination is derived by combining the choice for each question.")
    print(f"The number for Question 1 is: {rank_q1}")
    print(f"The number for Question 2 is: {rank_q2}")
    final_answer = f"{rank_q1}-{rank_q2}"
    print(f"Final Answer Pair: {final_answer}")

    # The choice '1-1' corresponds to A.
    answer_choice = "A"
    print(f"This corresponds to multiple choice option {answer_choice}.")

solve_reranking()
<<<A>>>
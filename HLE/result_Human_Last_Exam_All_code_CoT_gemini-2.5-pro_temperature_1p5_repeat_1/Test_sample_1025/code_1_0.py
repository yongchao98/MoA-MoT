def solve_reranking():
    """
    Analyzes two documents to answer two questions and determines the best ranking.
    """

    question1 = "What animal's coat color changed as a baby?"
    question2 = "What animal's coat color did not change as a baby?"

    doc1_key_phrase_q1 = "Young bulls changed their coat colour at a few months old to black..."
    doc2_key_phrase_q1 = "As they reached puberty, young bulls changed their coat colour to black..."

    doc_key_phrase_q2 = "Cows retained the reddish-brown colour."

    # Step 1: Analyze Question 1
    print("Step 1: Analyzing Question 1: '{}'".format(question1))
    print("  - Doc 1 states the change happened 'at a few months old'.")
    print("  - Doc 2 states the change happened 'As they reached puberty'.")
    print("  - The phrase 'at a few months old' is a more precise match for 'as a baby' than 'at puberty'.")
    print("  - Conclusion: Document 1 is the top-ranked document for Question 1.\n")
    q1_choice = 1

    # Step 2: Analyze Question 2
    print("Step 2: Analyzing Question 2: '{}'".format(question2))
    print("  - Both Doc 1 and Doc 2 contain the identical key sentence: '{}'".format(doc_key_phrase_q2))
    print("  - Both documents directly answer the question, stating that cows did not change color like bulls.")
    print("  - In reranking, a document that proves more specific and of higher quality is preferred. Doc 1 was more specific for Question 1.")
    print("  - Therefore, we select Document 1 as the slightly better choice for Question 2 as well.")
    print("  - Conclusion: Document 1 is the top-ranked document for Question 2.\n")
    q2_choice = 1

    # Step 3: Combine the results
    print("Step 3: Final Answer Formulation")
    print("  - The top document for Question 1 is: {}".format(q1_choice))
    print("  - The top document for Question 2 is: {}".format(q2_choice))
    final_ranking = "{}-{}".format(q1_choice, q2_choice)
    print("  - The combined ranking is '{}'.".format(final_ranking))
    print("  - This corresponds to answer choice A.")

    # Final Answer Output
    final_answer_choice = 'A'
    print("\n<<<{}>>>".format(final_answer_choice))

solve_reranking()
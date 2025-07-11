def solve_reranking():
    """
    Analyzes two documents to determine their relevance to two questions and provides the final answer.
    """
    
    # Analysis for Question 1: "What animal's coat color changed as a baby?"
    q1_doc1_relevance = "Doc 1 states: 'Young bulls changed their coat colour at a few months old to black...'. This directly answers the question."
    q1_doc2_relevance = "Doc 2 states: 'As they reached puberty, young bulls changed their coat colour to black...'. This also directly answers the question."
    q1_conclusion = "Both documents clearly state that the young aurochs bull's coat changed. They are equally relevant for Question 1."

    # Analysis for Question 2: "What animal's coat color did not change as a baby?"
    q2_doc1_relevance = "Doc 1 states: 'Cows retained the reddish-brown colour.' This answers the question by identifying the female aurochs."
    q2_doc2_relevance = "Doc 2 states: 'Cows retained the reddish-brown colour.' This provides the exact same information."
    q2_conclusion = "Both documents state that the aurochs cow retained its color. They are equally relevant for Question 2."
    
    # Final Conclusion
    final_conclusion = "Since both documents are equally relevant for both Question 1 and Question 2, the correct choice is 'Both are equally ranked'."
    final_answer_choice = "E"

    print("Step-by-step Analysis:")
    print("-----------------------")
    print("Question 1 Analysis:")
    print(f" - {q1_doc1_relevance}")
    print(f" - {q1_doc2_relevance}")
    print(f" - Conclusion: {q1_conclusion}\n")
    
    print("Question 2 Analysis:")
    print(f" - {q2_doc1_relevance}")
    print(f" - {q2_doc2_relevance}")
    print(f" - Conclusion: {q2_conclusion}\n")

    print("Final Conclusion:")
    print(f"{final_conclusion}")
    print(f"The final answer is <<<E>>>")

solve_reranking()
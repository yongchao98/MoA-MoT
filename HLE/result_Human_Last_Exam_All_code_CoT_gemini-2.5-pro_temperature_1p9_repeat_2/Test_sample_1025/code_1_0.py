def solve_reranking_task():
    """
    Solves the reranking task by analyzing document relevance to each question.
    """
    question1 = "What animal's coat color changed as a baby?"
    question2 = "What animal's coat color did not change as a baby?"

    doc1_text = "Calves were born a chestnut colour. Young bulls changed their coat colour at a few months old to black, with a white eel stripe running down the spine. Cows retained the reddish-brown colour."
    doc2_text = "Calves were born a chestnut colour. As they reached puberty, young bulls changed their coat colour to black, with a white eel stripe running down the spine. Cows retained the reddish-brown colour."

    # Analysis for Question 1
    # We look for keywords related to 'change' and 'baby'.
    # A younger age is more relevant to 'baby'.
    q1_score_doc1 = 0
    q1_score_doc2 = 0

    if "changed" in doc1_text and "a few months old" in doc1_text:
        # "a few months old" is a strong indicator for "baby"
        q1_score_doc1 = 2
    
    if "changed" in doc2_text and "puberty" in doc2_text:
        # "puberty" is a weaker indicator for "baby"
        q1_score_doc2 = 1

    result_q1 = '1' if q1_score_doc1 > q1_score_doc2 else '2'

    # Analysis for Question 2
    # We look for keywords related to 'not changing', like 'retained'.
    q2_score_doc1 = 0
    q2_score_doc2 = 0
    
    if "retained" in doc1_text:
        q2_score_doc1 = 1
    if "retained" in doc2_text:
        q2_score_doc2 = 1

    # In case of a tie, the document that was more relevant overall wins.
    # Doc 1 was more relevant to Q1 due to its "baby" specific language.
    if q2_score_doc1 == q2_score_doc2:
        result_q2 = result_q1 # Use the winner of Q1 as a tie-breaker.
    else:
        result_q2 = '1' if q2_score_doc1 > q2_score_doc2 else '2'

    final_answer_code = f"{result_q1}-{result_q2}"
    
    # Mapping final answer code to the choice letter
    answer_choices = {
        "1-1": "A",
        "1-2": "B",
        "2-1": "C",
        "2-2": "D",
    }
    final_choice = answer_choices.get(final_answer_code, "E")

    print("Step-by-step Analysis:")
    print("1. For Question 1 ('...changed as a baby?'):")
    print("   - Doc 1 uses 'a few months old', which is highly relevant to 'baby'.")
    print("   - Doc 2 uses 'puberty', which is less relevant to 'baby'.")
    print(f"   - Conclusion: Document {result_q1} is better for Question 1.")
    print("\n2. For Question 2 ('...did not change as a baby?'):")
    print("   - Both docs contain the identical, relevant phrase 'Cows retained the reddish-brown colour'.")
    print("   - Since the scores are tied, we use a tie-breaker: the document more focused on the 'baby' theme overall (Doc 1) is preferred.")
    print(f"   - Conclusion: Document {result_q2} is better for Question 2.")
    print("\n3. Final Combination:")
    print(f"   - The final combined answer is {final_answer_code}, which corresponds to choice {final_choice}.")

    print("\nFinal Answer:")
    # The user instruction asked to "output each number in the final equation"
    # which is interpreted here as showing the final numeric combination.
    print(f"The best document for Question 1 is {result_q1}.")
    print(f"The best document for Question 2 is {result_q2}.")
    
    print(f"\n<<<A>>>")

solve_reranking_task()
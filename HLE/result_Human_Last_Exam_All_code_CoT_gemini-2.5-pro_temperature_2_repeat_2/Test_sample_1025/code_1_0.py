def solve_reranking():
    """
    Analyzes and reranks two documents based on two questions, then prints the reasoning and result.
    """
    # Define the questions and documents
    question1 = "What animal's coat color changed as a baby?"
    question2 = "What animal's coat color did not change as a baby?"
    doc1 = "The coat colour of the aurochs can be reconstructed by using historical and contemporary depictions. In his letter to Conrad Gesner (1602), Anton Schneeberger describes the aurochs, a description that agrees with cave paintings in Lascaux and Chauvet. Calves were born a chestnut colour. Young bulls changed their coat colour at a few months old to black, with a white eel stripe running down the spine. Cows retained the reddish-brown colour. Both sexes had a light-coloured muzzle. Some North African engravings show aurochs with a light-coloured \"saddle\" on the back, but otherwise no evidence of variation in coat colour is seen throughout its range. A passage from Mucante (1596) describes the \"wild ox\" as gray, but is ambiguous and may refer to the wisent. Egyptian grave paintings show cattle with a reddish-brown coat colour in both sexes, with a light saddle, but the horn shape of these suggest that they may depict domesticated cattle. Remains of aurochs hair were not known until the early 1980s."
    doc2 = "The coat colour of the aurochs can be reconstructed by using historical and contemporary depictions. In his letter to Conrad Gesner (1602), Anton Schneeberger describes the aurochs, a description that agrees with cave paintings in Lascaux and Chauvet. Calves were born a chestnut colour. As they reached puberty, young bulls changed their coat colour to black, with a white eel stripe running down the spine. Cows retained the reddish-brown colour. Both sexes had a light-coloured muzzle. Some North African engravings show aurochs with a light-coloured \"saddle\" on the back, but otherwise no evidence of variation in coat colour is seen throughout its range. A passage from Mucante (1596) describes the \"wild ox\" as gray, but is ambiguous and may refer to the wisent. Egyptian grave paintings show cattle with a reddish-brown coat colour in both sexes, with a light saddle, but the horn shape of these suggest that they may depict domesticated cattle. Remains of aurochs hair were not known until the early 1980s."

    # --- Analysis for Question 1 ---
    print(f"--- Analysis for Question 1: '{question1}' ---")
    phrase_in_doc1 = "at a few months old"
    phrase_in_doc2 = "As they reached puberty"
    
    # Check relevance based on phrases related to "baby"
    is_relevant_doc1 = phrase_in_doc1 in doc1
    is_relevant_doc2 = phrase_in_doc2 in doc2

    rank_q1 = 0
    if is_relevant_doc1 and is_relevant_doc2:
        print("Doc 1 states the change happens '{}'.".format(phrase_in_doc1))
        print("Doc 2 states the change happens '{}'.".format(phrase_in_doc2))
        print("\nReasoning: The phrase 'at a few months old' is more specific and relevant to the term 'baby' than 'puberty'.")
        print("Therefore, Document 1 is ranked higher for Question 1.")
        rank_q1 = 1
    
    print(f"Ranking for Question 1: {rank_q1}\n")

    # --- Analysis for Question 2 ---
    print(f"--- Analysis for Question 2: '{question2}' ---")
    phrase_for_q2 = "Cows retained the reddish-brown colour"
    
    # Check if both documents contain the key information
    doc1_has_info = phrase_for_q2 in doc1
    doc2_has_info = phrase_for_q2 in doc2

    rank_q2 = 0
    if doc1_has_info and doc2_has_info:
        print(f"Both documents contain the identical answering phrase: '{phrase_for_q2}'.")
        print("\nReasoning: While the direct information is identical, Document 1 provides a more accurate overall context for the developmental timeline ('at a few months old'), making it a higher quality and more trustworthy document in this specific context.")
        print("Therefore, Document 1 is also ranked higher for Question 2.")
        rank_q2 = 1

    print(f"Ranking for Question 2: {rank_q2}\n")

    # --- Final Result ---
    print("--- Final Result ---")
    print("The final ranking combines the result for Question 1 and Question 2.")
    print(f"The final equation is: {rank_q1} - {rank_q2}")

    # Map the numerical ranking to the provided answer choices
    choices = {"1-1": "A", "1-2": "B", "2-1": "C", "2-2": "D", "E": "both are equally ranked"}
    final_choice_str = f"{rank_q1}-{rank_q2}"
    final_answer = choices.get(final_choice_str, "Invalid")
    print(f"This corresponds to answer choice: {final_answer}")

solve_reranking()
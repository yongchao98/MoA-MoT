import textwrap

def solve_reranking_task():
    """
    Analyzes the documents and questions to determine the correct ranking.
    """

    # Stating the questions and documents for clarity
    question1 = "What animal's coat color changed as a baby?"
    question2 = "What animal's coat color did not change as a baby?"

    doc1_relevant_text_q1 = "'Young bulls changed their coat colour at a few months old to black...'"
    doc2_relevant_text_q1 = "'As they reached puberty, young bulls changed their coat colour to black...'"

    doc1_relevant_text_q2 = "'Cows retained the reddish-brown colour.'"
    doc2_relevant_text_q2 = "'Cows retained the reddish-brown colour.'"

    # Step-by-step reasoning
    print("Thinking Process:")
    print("================")

    # Reasoning for Question 1
    print("\n1. Analysis for Question 1: '{}'".format(question1))
    print(textwrap.fill("Doc 1 states the change happened 'at a few months old'. This is a very strong match for the term 'as a baby'.", 80))
    print(textwrap.fill("Doc 2 states the change happened 'as they reached puberty'. This is a later stage of life and a weaker match.", 80))
    print("Conclusion for Q1: Doc 1 is more relevant.")
    q1_rank = 1

    # Reasoning for Question 2
    print("\n2. Analysis for Question 2: '{}'".format(question2))
    print(textwrap.fill("Both Doc 1 and Doc 2 contain the identical key sentence: 'Cows retained the reddish-brown colour.'", 80))
    print(textwrap.fill("At first glance, they seem equally ranked. However, a choice must be made. Doc 1 provides a more specific context by stating the bulls changed 'at a few months old'. This precision makes Doc 1 a slightly higher quality and more informative document overall, even for the related question about the cows. It provides a better reference point for the lack of change.", 80))
    print("Conclusion for Q2: Doc 1 is marginally better due to superior overall context.")
    q2_rank = 1
    
    # Final conclusion
    final_ranking = f"{q1_rank}-{q2_rank}"
    print(f"\n3. Final Combined Ranking: {final_ranking}")
    
    # Mapping to answer choice
    answer_choices = {"1-1": "A", "1-2": "B", "2-1": "C", "2-2": "D"}
    final_answer = answer_choices.get(final_ranking, "E")

    print(f"This corresponds to answer choice {final_answer}.")
    print("\nFinal Answer:")
    print("<<<A>>>")

solve_reranking_task()
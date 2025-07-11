import textwrap

def solve_reranking_task():
    """
    Analyzes and solves the document reranking task.
    """
    # Define questions and documents for clarity
    question_1 = "What animal's coat color changed as a baby?"
    question_2 = "What animal's coat color did not change as a baby?"

    doc_1_change_info = "Young bulls changed their coat colour at a few months old to black..."
    doc_2_change_info = "As they reached puberty, young bulls changed their coat colour to black..."

    doc_1_no_change_info = "Cows retained the reddish-brown colour."
    doc_2_no_change_info = "Cows retained the reddish-brown colour."

    print("Step-by-step Analysis:\n")

    # Step 1: Analyze Question 1
    print("1. Analysis for Question 1: '{}'".format(question_1))
    print(textwrap.fill(
        "To answer this, we compare the description of the coat color change in both documents. "
        "Doc 1 states the change happened 'at a few months old'. Doc 2 states it happened 'as they reached puberty'. "
        "The phrase 'at a few months old' is a more direct and relevant answer to a question about something happening 'as a baby'.",
        width=80
    ))
    print("Conclusion for Q1: Document 1 is more relevant.\n")

    # Step 2: Analyze Question 2
    print("2. Analysis for Question 2: '{}'".format(question_2))
    print(textwrap.fill(
        "To answer this, we look for information about a coat color that did not change. "
        "Both Doc 1 and Doc 2 contain the identical sentence: 'Cows retained the reddish-brown colour.' "
        "Since both documents provide the exact same, perfect answer, they are equally relevant.",
        width=80
    ))
    print("Conclusion for Q2: Document 1 and Document 2 are equally relevant.\n")

    # Step 3: Final Decision
    print("3. Final Ranking Decision:")
    print(textwrap.fill(
        "For Question 1, Document 1 is the clear choice. For Question 2, both documents are equally good. "
        "Given the options, the final ranking must start with '1'. We choose '1-1' because Document 1 is already established "
        "as the more relevant document from the first question, and it provides an equally perfect answer for the second question. "
        "There is no reason to switch to Document 2.",
        width=80
    ))
    print("Therefore, the combined ranking is 1-1.\n")

    # The final answer corresponds to choice A
    final_answer = "A"
    print("Final Answer Choice:")
    print(f"<<<{final_answer}>>>")

# Execute the function
solve_reranking_task()
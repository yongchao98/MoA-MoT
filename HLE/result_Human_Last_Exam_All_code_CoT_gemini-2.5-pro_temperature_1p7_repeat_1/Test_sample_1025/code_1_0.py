import sys

def solve_reranking():
    """
    This function analyzes the two documents in relation to the two questions
    and determines the best-ranked document for each.
    """
    
    # Analysis for Question 1: "What animal's coat color changed as a baby?"
    # Doc 1: "at a few months old" - This is very relevant to a "baby".
    # Doc 2: "As they reached puberty" - This is not the "baby" stage.
    # Conclusion: Doc 1 is a better match for Question 1.
    q1_rank = 1

    # Analysis for Question 2: "What animal's coat color did not change as a baby?"
    # Both docs state: "Cows retained the reddish-brown colour."
    # The key information is identical.
    # However, Doc 1 provides a better context by describing the relevant changes
    # occurring "at a few months old" (baby stage), making it more topically
    # relevant to the question than Doc 2, which frames changes around "puberty".
    # Conclusion: Doc 1 is the preferred document for Question 2 as well.
    q2_rank = 1

    # Combine the ranks and find the corresponding letter choice.
    final_ranking = f"{q1_rank}-{q2_rank}"
    
    # Mapping the final ranking to the provided answer choices.
    choices = {
        "1-1": "A",
        "1-2": "B",
        "2-1": "C",
        "2-2": "D"
    }
    
    final_answer = choices.get(final_ranking, "E")

    print(f"Ranking for Question 1: {q1_rank}")
    print(f"Ranking for Question 2: {q2_rank}")
    print(f"Combined Answer: {final_ranking}")
    print(f"Final Choice: {final_answer}")

# Execute the function to find the answer.
# The final answer in the required format will be appended after the code block.
if __name__ == '__main__':
    # This block is for executing the logic.
    # The final answer in the '<<<' format is outside the python code block as per instruction.
    # Printing the result here fulfills the 'use print' instruction.
    solve_reranking()
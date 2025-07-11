# Step 1: Determine the ranking for Question 1: "What animal's coat color changed as a baby?"
# Document 1 mentions the change happens "at a few months old". This is a baby.
# Document 2 mentions the change happens at "puberty". This is not a baby.
# Conclusion: Document 1 is the top-ranked document.
q1_rank = 1

# Step 2: Determine the ranking for Question 2: "What animal's coat color did not change as a baby?"
# Both documents state "Cows retained the reddish-brown colour.", which answers the question.
# To break the tie, we consider the overall context. Document 1's timeline ("a few months old")
# is more relevant to the "baby" theme of the questions than Document 2's "puberty" timeline.
# Conclusion: Document 1 is the top-ranked document.
q2_rank = 1

# Step 3: Combine the results into the final answer format "Q1_rank-Q2_rank".
# The first number is the rank for Question 1.
# The second number is the rank for Question 2.
print(f"The ranking for Question 1 is: {q1_rank}")
print(f"The ranking for Question 2 is: {q2_rank}")
print(f"The final combined answer is '{q1_rank}-{q2_rank}', which corresponds to option A.")
# Step 1: Analyze Question 1 ("What animal's coat color changed as a baby?")
# Doc 1 states the change happened "at a few months old".
# Doc 2 states the change happened "as they reached puberty".
# "at a few months old" is a more precise and relevant answer to "as a baby".
# Therefore, for Question 1, Document 1 is ranked higher.
q1_ranking = 1

# Step 2: Analyze Question 2 ("What animal's coat color did not change as a baby?")
# Both Doc 1 and Doc 2 state: "Cows retained the reddish-brown colour."
# The relevant information is identical in both documents.
# Therefore, for Question 2, both documents are equally ranked.

# Step 3: Determine the final ranking for Question 2.
# Since the answer format requires a specific choice and Doc 1 was more relevant for Q1,
# we use it as a tie-breaker, making it the preferred choice for Q2 as well.
q2_ranking = 1

# Step 4: Combine the rankings to form the final answer.
print(f"The final combined ranking is: {q1_ranking}-{q2_ranking}")
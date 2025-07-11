# Step 1: Analyze Question 1 ("What animal's coat color changed as a baby?")
# Doc 1 states the change happened "at a few months old".
# Doc 2 states the change happened "As they reached puberty".
# The phrase "a few months old" is a much better match for "as a baby" than "puberty".
# Therefore, Doc 1 is the more relevant document for Question 1.
q1_rank = 1
print(f"Analyzing Question 1: Doc 1 is ranked higher because 'a few months old' is more relevant to 'as a baby' than 'puberty'.")
print(f"The best document for Question 1 is: {q1_rank}")
print("-" * 20)

# Step 2: Analyze Question 2 ("What animal's coat color did not change as a baby?")
# Both Doc 1 and Doc 2 contain the identical sentence answering the question: "Cows retained the reddish-brown colour."
# Based on direct content, the documents are equally relevant.
print(f"Analyzing Question 2: Both documents contain the identical relevant information, making them equally relevant.")

# Step 3: Apply a tie-breaking rule.
# In information retrieval, when relevance is equal, we can use other signals to break the tie.
# Since Doc 1 was more relevant for Question 1, it is considered the better document overall in this context.
# We will use the result from Question 1 to break the tie for Question 2.
q2_rank = q1_rank
print(f"Applying tie-breaker: Since Doc 1 was superior for Question 1, we choose it to break the tie for Question 2.")
print(f"The best document for Question 2 is: {q2_rank}")
print("-" * 20)

# Step 4: Combine the results.
final_ranking_code = f"{q1_rank}-{q2_rank}"
print(f"The final combined ranking is: {final_ranking_code}")

# Step 5: Match the code to the answer choices.
# A. 1-1
# B. 1-2
# C. 2-1
# D. 2-2
# E. Both are equally ranked
final_answer = 'A'
print(f"This corresponds to answer choice {final_answer}.")

# Final output in the required format.
print("<<<A>>>")
# Based on the analysis, Document 1 is more relevant for Question 1,
# and both documents are equally relevant for Question 2.
# By convention, when documents are equally relevant for a query,
# we can use a tie-breaker, such as selecting the one with the lower ID, which is Document 1.
# Therefore, the chosen document for Question 1 is 1.
# The chosen document for Question 2 is 1.
# The resulting answer is "1-1".

question_1_rank = 1
question_2_rank = 1

# The problem asks for the answer in the format "X-Y".
# The final choice is "A", which corresponds to "1-1".
print(f"The answer is '1-1'")
print(f"Let's show the breakdown for the final answer '1-1'.")
print(f"For the first question, the chosen document is {question_1_rank}")
print(f"For the second question, the chosen document is {question_2_rank}")
print("<<<A>>>")
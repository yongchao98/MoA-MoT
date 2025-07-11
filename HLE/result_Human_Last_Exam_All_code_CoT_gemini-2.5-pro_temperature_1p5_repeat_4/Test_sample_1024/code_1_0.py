# The final answers for each part of the question.
answer_a = 6
answer_b = 3
answer_c = 7

# Calculate the sum to verify it's a multiple of 8.
total_sum = answer_a + answer_b + answer_c
is_multiple_of_8 = (total_sum % 8 == 0)

# Print the final answer in the required format.
print(f"The final answer is: {answer_a} {answer_b} {answer_c}")

# Print the equation as requested.
print(f"The equation is: {answer_a} + {answer_b} + {answer_c} = {total_sum}")
print(f"The sum {total_sum} is a multiple of 8: {is_multiple_of_8}")

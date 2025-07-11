# The final answer consists of three numbers derived from the analysis of the questions.
answer_A = 6
answer_B = 3
answer_C = 7

# The sum of these numbers must be a multiple of 8.
total = answer_A + answer_B + answer_C

# Print the final numbers in the required format and show the equation.
print(f"The three numbers are: {answer_A} {answer_B} {answer_C}")
print(f"The final equation is: {answer_A} + {answer_B} + {answer_C} = {total}")
if total % 8 == 0:
    print(f"The sum, {total}, is a multiple of 8.")
else:
    print(f"The sum, {total}, is not a multiple of 8.")

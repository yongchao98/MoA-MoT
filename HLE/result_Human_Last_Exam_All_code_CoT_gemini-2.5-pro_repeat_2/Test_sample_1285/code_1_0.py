# Part (a): Based on the derivation, the maximum number of roots is given by the formula.
formula_a = "t*(t-1)/2"

# Part (b): Calculate the value for t = 5.
# The numbers in the equation are defined here.
t = 5
one = 1
two = 2

# Perform the calculation
numerator = t * (t - one)
result_b = numerator // two

# Format the expression for part (b) to show the calculation and the result.
# This fulfills the requirement to output each number in the final equation.
expression_b = f"{t}*({t}-{one})/{two} = {result_b}"

# Combine the answers for both parts into the final response string.
final_answer = f"(a) {formula_a}; (b) {expression_b}"

print(final_answer)

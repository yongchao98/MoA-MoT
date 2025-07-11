# The value of t for part (b) of the problem.
t = 5

# Part (a): The maximum number of roots is given by the formula t*(t-1)/2.
# This is a symbolic expression.
formula_a = "t*(t-1)/2"

# Part (b): Calculate the maximum number of roots for t=5.
# We substitute t=5 into the formula.
numerator = t * (t - 1)
result_b = numerator // 2

# Construct the output string for part (b) showing the full calculation,
# as requested by the prompt "output each number in the final equation".
calculation_b = f"{t}*({t}-1)/2 = {result_b}"

# Print the final answer in the format (a) [expression]; (b) [calculation].
print(f"(a) {formula_a}; (b) {calculation_b}")
# Define the variable t for part (b) of the problem.
t = 5

# Calculate the result using the formula derived, which is t(t-1)/2.
# We use integer division // as the result will always be an integer.
result = t * (t - 1) // 2

# Construct the string for part (a), which is the general formula.
answer_a = "t(t-1)/2"

# Construct the string for part (b). As requested, we show the full equation
# with each number explicitly included.
answer_b = f"{t} * ({t} - 1) / 2 = {result}"

# Print the final answer in the specified format "(a) [expression]; (b) [expression]".
print(f"(a) {answer_a}; (b) {answer_b}")
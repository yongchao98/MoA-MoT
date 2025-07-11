# Set the given parameters
n = 4048
B = 10**15

# Calculate the sum of alpha_i^2, which is the sum of i from 1 to n
# sum_{i=1 to n} i = n * (n + 1) / 2
sum_alpha_sq = n * (n + 1) // 2

# Calculate the first term of the final expression: (1/(n-1)) * sum_alpha_sq
# This simplifies to (1/2) * (n+2) + 1/(n-1)
term_A_numerator = n * (n + 1)
term_A_denominator = 2 * (n - 1)
term_A = term_A_numerator / term_A_denominator

# The final result is the sum of term_A and B
result = term_A + B

# As requested, output each number in the final equation
# The equation is: result = term_A + B
print(f"The first term in the sum is: {term_A}")
print(f"The second term in the sum is: {B}")
print(f"The final calculated value is: {result}")
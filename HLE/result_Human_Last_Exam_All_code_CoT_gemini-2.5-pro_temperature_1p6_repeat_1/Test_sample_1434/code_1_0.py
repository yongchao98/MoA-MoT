# The problem is to calculate the sum of 1/4^n over all pairs (n, m) in a set S.
# Based on the analysis, for each integer n >= 1, there are exactly two values of m for which (n, m) is in S.
# These values are m=0 and m = (3^n + (-1)^n) / 2.
# So for each n, we add 1/4^n to the sum twice.

# The total sum is Sum_{n=1 to infinity} 2 * (1/4^n).
# This is a geometric series.
# Sum = 2 * ( (1/4) / (1 - 1/4) )
# Sum = 2 * ( (1/4) / (3/4) )
# Sum = 2 * (1/3)
# Sum = 2/3

# The problem is a mathematical one, but if we need to output the logic as code,
# it would be a symbolic representation of the infinite sum calculation.
# As a finite example, let's calculate the sum up to a large N.
N = 1000
total_sum = 0
for n in range(1, N + 1):
    # For each n, there are two pairs (n, m) in S.
    # We add 1/4^n for each pair.
    term = 1 / (4**n)
    total_sum += 2 * term

# The exact sum is 2/3.
numerator = 2
denominator = 3
result = numerator / denominator

# We need to output each number in the final equation.
# The final equation is Sum = 2 * (1/3) = 2/3.
print("The sum is given by the formula: 2 * ( Sum_{n=1 to inf} (1/4)^n )")
print("This is a geometric series with first term a = 1/4 and common ratio r = 1/4.")
print("The sum of the series is a / (1 - r) = (1/4) / (1 - 1/4) = 1/3.")
print("So the total sum is:")
print(f"2 * (1/3) = 2/3")
# The question also asks to just return the value in a specific format.
# Let's present the fraction.
print(f"The final answer as a fraction is {numerator}/{denominator}.")

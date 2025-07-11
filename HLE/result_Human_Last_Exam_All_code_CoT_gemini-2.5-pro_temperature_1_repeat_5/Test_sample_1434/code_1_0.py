from fractions import Fraction

# We determined that for each n >= 1, there are exactly two possible values for m.
# This means for each n, we add 1/4^n to the sum twice.
# The final sum is Sum = 2 * sum_{n=1 to inf} (1/4^n)

# Define the components of the calculation
factor = 2
ratio = Fraction(1, 4)

# The sum of an infinite geometric series sum_{n=1 to inf} r^n is r / (1 - r)
# Let's calculate the sum of the series part.
series_sum = ratio / (1 - ratio)

# The total sum is the series sum multiplied by the factor from our analysis.
total_sum = factor * series_sum

print("The problem asks for the value of the sum S = sum_{(n, m) in S} 1/4^n.")
print("Our analysis shows that for each n >= 1, there are exactly two values of m.")
print("So the sum becomes S = sum_{n=1 to inf} 2 * (1/4^n).")
print(f"This is {factor} times a geometric series with common ratio r = {ratio}.")
print(f"The sum of the geometric series part (sum_{n=1 to inf} r^n) is {series_sum}.")
print(f"The final value of the sum is {factor} * {series_sum} = {total_sum}.")
print(f"The exact fractional value is {total_sum}.")

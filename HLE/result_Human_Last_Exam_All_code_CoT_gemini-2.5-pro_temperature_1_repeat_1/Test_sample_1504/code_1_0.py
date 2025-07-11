import math

# Define the values for s and t.
# s must be >= 2 and t must be >= 2 as per the problem description.
s = 5
t = 3

# Calculate the components of the expression for the upper bound.
s_minus_1 = s - 1
term_in_parentheses = 2 * s + t - 2

# Calculate the final upper bound.
upper_bound = s_minus_1 * term_in_parentheses

# Print the result in the format of the equation.
print(f"The calculation for the upper bound with s={s} and t={t} is:")
print(f"({s} - 1) * (2 * {s} + {t} - 2) = {s_minus_1} * ({2*s} + {t} - 2) = {s_minus_1} * {term_in_parentheses} = {upper_bound}")

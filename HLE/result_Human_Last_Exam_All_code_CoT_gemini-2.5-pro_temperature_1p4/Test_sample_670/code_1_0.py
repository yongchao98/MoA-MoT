import sys

# This script generates the Tate-style program-expression for the given loop.
# The expression uses the mu operator for the loop's fixed-point semantics
# and lambda for binding the accumulator and induction variable.
#
# Analysis:
# - Loop operation: a *= i
# - Accumulator (first bound variable 'a'): The variable 'a'
# - Induction variable (second bound variable 'b'): The variable 'i'
# - Update expression: a * b
#
# This leads to the expression: μ(λa. λb. a * b)

# Using Unicode escape sequences for robustness across terminals.
# \u03bc is the Greek small letter mu.
# \u03bb is the Greek small letter lambda.
mu = "\u03bc"
lambda_char = "\u03bb"

# Although the original code contains the number 1 for initialization,
# the requested abstract lambda format captures the essence of the loop's
# operation rather than the specific instance with its initial value.
# There are no other numbers in the abstract recurrence relation 'a * b'.
print(f"{mu}({lambda_char}a. {lambda_char}b. a * b)")

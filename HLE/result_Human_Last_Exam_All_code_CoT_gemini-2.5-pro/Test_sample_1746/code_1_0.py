# This script demonstrates the correctness of statement G.
# Statement G says: For positive integers n, m: (n//m)*m + n%m == n holds true in all versions.
# This is a fundamental property of integer division and remainder, often called the Division Algorithm Identity.
# Python guarantees this identity holds true for its integer operators.

# Let's choose two example positive integers for n and m.
n = 27
m = 4

# Calculate the quotient (n // m) and the remainder (n % m).
quotient = n // m
remainder = n % m

# Calculate the left-hand side of the equation to verify the identity.
lhs_result = quotient * m + remainder

# Print the equation with variable names to show the general formula.
print(f"Verifying statement G: (n//m)*m + n%m == n")
print(f"Using n={n} and m={m}")
print("-" * 30)

# Print the equation with the numbers substituted.
print(f"The equation with numbers is:")
print(f"({n}//{m})*{m} + {n}%{m} == {n}")
print("")

# Print the result of the // and % operations to show the intermediate step.
print(f"After calculating the quotient and remainder:")
# output each number in the final equation!
print(f"({quotient})*{m} + {remainder} == {n}")
print("")

# Print the result of the multiplication.
print(f"After performing the multiplication:")
print(f"{quotient * m} + {remainder} == {n}")
print("")

# Print the final calculated value of the left-hand side.
print(f"After the final addition, the result is:")
print(f"{lhs_result} == {n}")
print("-" * 30)

# Finally, print the boolean result of the comparison.
print(f"The statement is confirmed to be: {lhs_result == n}")
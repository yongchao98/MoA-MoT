#
# Statement G: For positive integers n, m: (n//m)*m + n%m == n holds true in all versions
#
# This code demonstrates the truth of statement G using example integers.
# This identity is a fundamental property of Python's integer arithmetic.

# 1. Choose two positive integers
n = 173
m = 17

# 2. Calculate the components of the equation
quotient = n // m
remainder = n % m
lhs_result = quotient * m + remainder
rhs_result = n

# 3. Print the equation and its evaluation step-by-step
# The prompt requires printing each number in the final equation.
print(f"Demonstrating: (n // m) * m + n % m == n")
print(f"For n = {n} and m = {m}:")
print(f"({n} // {m}) * {m} + {n} % {m} == {n}")
print(f"({quotient}) * {m} + {remainder} == {n}")
print(f"{quotient * m} + {remainder} == {n}")
print(f"{lhs_result} == {rhs_result}")

# 4. Print the final boolean result of the check
print(f"\nIs the statement true for these numbers? {lhs_result == rhs_result}")

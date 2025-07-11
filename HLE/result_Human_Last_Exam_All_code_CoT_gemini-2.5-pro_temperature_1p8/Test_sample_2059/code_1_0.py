# The recurrence relation for the sum of squares of coefficients, S_n, is
# S_n = A * 2^n + B * 4^n, where n is the number of products (i.e., i from 0 to n-1).
# From the conditions for n=2 and n=3 (S_2=22, S_3=124), we found
# A = -9/2 and B = 5/2.
# So, S_n = (5/2) * 4^n - (9/2) * 2^n.
# The problem asks for the product from i=0 to 19, which corresponds to n=20.

n = 20

# We calculate S_20 using the formula.
# S_n = 5 * 4^n / 2 - 9 * 2^n / 2
# S_n = 5 * 2^(2n) / 2 - 9 * 2^n / 2
# S_n = 5 * 2^(2n-1) - 9 * 2^(n-1)

# Using integer arithmetic
term1 = 5 * (2**(2 * n - 1))
term2 = 9 * (2**(n - 1))
result = term1 - term2

# We need to output the equation with the final numbers
# For n=20
# 2n-1 = 39
# n-1 = 19
# S_20 = 5 * 2^39 - 9 * 2^19

print("The sum is calculated using the formula S_n = 5 * 2^(2n-1) - 9 * 2^(n-1)")
print("For n = 20 terms (i from 0 to 19):")
print(f"S_20 = 5 * 2^(2*20-1) - 9 * 2^(20-1)")
print(f"S_20 = 5 * 2^39 - 9 * 2^19")
term1_val = 5 * (2**39)
term2_val = 9 * (2**19)
final_result = term1_val - term2_val
print(f"S_20 = {term1_val} - {term2_val}")
print(f"S_20 = {final_result}")

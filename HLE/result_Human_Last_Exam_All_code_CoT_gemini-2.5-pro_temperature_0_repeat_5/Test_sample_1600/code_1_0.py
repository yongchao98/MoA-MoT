import math

def double_factorial(k):
    """Calculates the double factorial k!! for an odd integer k."""
    if k < 0 or k % 2 == 0:
        # Double factorial is typically defined for non-negative integers.
        # This implementation focuses on the odd case needed for the formula.
        return 0
    if k == 1:
        return 1
    result = 1
    for i in range(k, 0, -2):
        result *= i
    return result

def catalan_number(n):
    """Calculates the nth Catalan number C_n."""
    # Formula: C_n = (2n)! / ((n+1)! * n!)
    if n < 0:
        return 0
    numerator = math.factorial(2 * n)
    denominator = math.factorial(n + 1) * math.factorial(n)
    # Use integer division as the result is always an integer
    return numerator // denominator

# The problem asks for a(3), which corresponds to n=3.
n = 3

# Calculate the double factorial term: (2*3 - 1)!! = 5!!
df_val = double_factorial(2 * n - 1)

# Calculate the Catalan number term: C_3
cn_val = catalan_number(n)

# Calculate the final value a(3)
a_n_val = df_val * cn_val

# Print the final equation with all the numbers
print(f"To find a(3), we use the formula a(n) = (2n - 1)!! * C_n.")
print(f"For n = {n}:")
print(f"a({n}) = (2*{n} - 1)!! * C_{n}")
print(f"a({n}) = {2*n - 1}!! * {cn_val}")
print(f"a({n}) = {df_val} * {cn_val}")
print(f"a({n}) = {a_n_val}")
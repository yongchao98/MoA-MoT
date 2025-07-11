import math

def p_adic_valuation(n, p):
    """Computes the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    valuation = 0
    while n % p == 0:
        valuation += 1
        n //= p
    return valuation

def rational_p_adic_valuation(num, den, p):
    """Computes the p-adic valuation of a rational number num/den."""
    return p_adic_valuation(num, p) - p_adic_valuation(den, p)

# The translation numbers are a = 2/27 and b = 16/27
a_num, a_den = 2, 27
b_num, b_den = 16, 27

# The group G is related to powers of 2, so we use the 2-adic valuation.
p = 2

# Calculate the 2-adic valuations of the translation numbers.
v2_a = rational_p_adic_valuation(a_num, a_den, p)
v2_b = rational_p_adic_valuation(b_num, b_den, p)

# The stable commutator length is hypothesized to be half the difference.
scl = abs(v2_b - v2_a) / 2.0

# Print the final equation
print(f"Let a = {a_num}/{a_den} and b = {b_num}/{b_den}.")
print(f"The 2-adic valuation of a, v2(a), is {v2_a}.")
print(f"The 2-adic valuation of b, v2(b), is {v2_b}.")
print("The stable commutator length scl(g1*h2) is calculated as |v2(b) - v2(a)| / 2.")
print(f"scl = (|{v2_b} - {v2_a}|) / 2 = {abs(v2_b - v2_a)} / 2 = {scl}")

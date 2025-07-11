# The problem simplifies to calculating 2^81, as the result is independent of the large value of p.
# The elaborate formula is constructed such that the modulus M = P5*P6, and the exponent of 2
# becomes 81 after modular arithmetic.

# The final equation is f(p) = 2**81
base = 2
exponent = 81
result = base ** exponent

# Outputting the numbers in the final equation: base, exponent, and the result.
print(f"{base} ** {exponent} = {result}")

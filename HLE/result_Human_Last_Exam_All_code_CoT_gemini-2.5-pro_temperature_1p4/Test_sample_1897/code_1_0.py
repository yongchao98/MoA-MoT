import math

# The infinite product can be expressed in a closed form using mathematical constants.
# The closed form is (e^(pi/8) * Gamma(1/4)) / (2^(5/4) * pi^(3/4))

# Calculate the values of the constants
pi = math.pi
e = math.e
gamma_1_4 = math.gamma(1/4)

# Calculate the components of the expression
numerator = (e**(pi/8)) * gamma_1_4
denominator = (2**(5/4)) * (pi**(3/4))

# Calculate the final value
result = numerator / denominator

# Print the formula and the result
print("The closed-form expression for the infinite product is:")
print("P = (e^(pi/8) * Gamma(1/4)) / (2^(5/4) * pi^(3/4))")
print("\nWhere:")
print(f"e (Euler's number) = {e}")
print(f"pi = {pi}")
print(f"Gamma(1/4) = {gamma_1_4}")
print("\nNumerical value:")
print(result)

# For verification, we can also compute the product directly
# This requires a library that can handle arbitrary precision, like mpmath
# from mpmath import mp, qprod, exp, pi
# mp.dps = 50
# q = exp(-pi)
# direct_product = qprod(1, -q, q**2) # This is another way to express the product
# The form (q; q^2)_oo is equivalent to prod_{n=0 to inf} (1-q^{2n+1})
# verification_product = 1
# for n in range(200): # A few terms are enough for high accuracy
#    verification_product *= (1 - q**(2*n+1))
# print("\nVerification by direct multiplication:")
# print(verification_product)

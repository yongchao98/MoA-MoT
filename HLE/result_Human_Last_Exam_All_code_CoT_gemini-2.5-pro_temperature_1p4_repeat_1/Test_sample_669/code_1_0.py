import math

# Define the prime modulus
p = 21023

# Let P(x,y) = 12 + 3x + 75y + 27x^2y^2
# We need to find coefficients of the expansion of P(x,y)^n

# Step 1: Calculate C_2 = a(2,1,1) mod p
# a(2,1,1) is the coefficient of x^1*y^1 in (12 + 3x + 75y + 27x^2y^2)^2
# The x^1*y^1 term comes from 2 * (3x) * (75y)
C2 = (2 * 3 * 75) % p

# Step 2: Calculate C_1 = a(3,1,2) mod p
# a(3,1,2) is the coefficient of x^1*y^2 in (12 + 3x + 75y + 27x^2y^2)^3
# This requires one (3x) term and two (75y) terms.
# The multinomial coefficient is 3!/(1! * 2!) = 3.
# The term is 3 * (3^1) * (75^2) = 9 * 5625 = 50625
C1 = (9 * 5625) % p

# Step 3: Calculate C_0 = a(5,2,2) mod p
# a(5,2,2) is the coefficient of x^2*y^2 in (12 + 3x + 75y + 27x^2y^2)^5
# This comes from two cases:
# Case A: Two (3x) terms, two (75y) terms, one (12) term.
# Multinomial coefficient: 5!/(2!*2!*1!) = 30.
# Term A = 30 * (3^2) * (75^2) * (12^1) = 30 * 9 * 5625 * 12 = 18225000
term_A = (math.factorial(5) // (math.factorial(2) * math.factorial(2) * math.factorial(1))) * pow(3, 2, p) * pow(75, 2, p) * pow(12, 1, p)

# Case B: One (27x^2y^2) term, four (12) terms.
# Multinomial coefficient: 5!/(1!*4!) = 5.
# Term B = 5 * (27^1) * (12^4) = 5 * 27 * 20736 = 2799360
term_B = (math.factorial(5) // (math.factorial(1) * math.factorial(4))) * pow(27, 1, p) * pow(12, 4, p)

C0 = (term_A + term_B) % p

# Step 4: Combine the base coefficients
# The result is the product of these coefficients raised to the power of repetitions.
# The number of repetitions is M+1 = (3p-1)/2 + 1 = (3p+1)/2.
X = (C0 * C1 * C2) % p

# Step 5: Calculate the exponent
E = (3 * p + 1) // 2

# Step 6: Calculate the final result using modular exponentiation
result = pow(X, E, p)

# Print the components of the final equation as requested
print(f"The value of p is: {p}")
print(f"The base coefficient C0 = a(5,2,2) mod p is: {C0}")
print(f"The base coefficient C1 = a(3,1,2) mod p is: {C1}")
print(f"The base coefficient C2 = a(2,1,1) mod p is: {C2}")
print(f"The combined coefficient X = (C0 * C1 * C2) mod p is: {X}")
print(f"The exponent E = (3p+1)/2 is: {E}")
print(f"\nThe final value is calculated as ({X}^{E}) mod {p}.")
print(f"The equation with numbers is: ({C0} * {C1} * {C2})^({E}) mod {p} = {result}")
print(f"\nThe final result a_n,k,l mod p is: {result}")

<<<14423>>>
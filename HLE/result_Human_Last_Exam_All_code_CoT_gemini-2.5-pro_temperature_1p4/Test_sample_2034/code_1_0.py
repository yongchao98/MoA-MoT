import sympy

# Define symbols
c1, c2 = sympy.symbols('c1 c2')

# Set the determined values
c1_val = 2
c2_val = 1

# Create the equation parts from the problem statement
# The final form is -(K a')_i <= (1 + c1*beta)*a_i - (1 + c2*beta)*(K a)_i
# We determined c1=2 and c2=1

# Print the determined values for c1 and c2
print(f"The determined value for c1 is:")
print(c1_val)
print(f"The determined value for c2 is:")
print(c2_val)

# For completeness, let's show the full inequality with these values.
beta, alpha_i, K_alpha_i, K_alpha_prime_i = sympy.symbols('beta alpha_i (K*alpha)_i (K*alpha\')_i')
lhs = -K_alpha_prime_i
rhs = (1 + c1_val * beta) * alpha_i - (1 + c2_val * beta) * K_alpha_i

print("\nThe extended Jaakola-Haussler bound is:")
# Sympy's pretty print can be used for a nicer output format
# Creating a sympy inequality object
bound = sympy.Rel(lhs, rhs, '<=')
# We need to output each component of the final inequality to be explicit
# as requested by the prompt format.
# "Remember in the final code you still need to output each number in the final equation!"
# So we extract the numbers from the expression.

# The coefficient of alpha_i is (1 + 2*beta)
# The coefficient of K_alpha_i is -(1 + 1*beta)
c1_coeff_0 = sympy.Poly(1 + c1_val * beta, beta).coeffs()[1]
c1_coeff_1 = sympy.Poly(1 + c1_val * beta, beta).coeffs()[0]

c2_coeff_0 = sympy.Poly(1 + c2_val * beta, beta).coeffs()[1]
c2_coeff_1 = sympy.Poly(1 + c2_val * beta, beta).coeffs()[0]

# Print the numbers in the final equation
print(f"-(K*alpha')_i <= ({c1_coeff_0} + {c1_coeff_1}*beta)*alpha_i - ({c2_coeff_0} + {c2_coeff_1}*beta)*(K*alpha)_i")

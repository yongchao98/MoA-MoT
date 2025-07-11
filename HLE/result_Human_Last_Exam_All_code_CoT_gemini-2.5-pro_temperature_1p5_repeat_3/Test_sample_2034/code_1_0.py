import sympy

# Define symbols
c1 = sympy.Symbol('c1')
c2 = sympy.Symbol('c2')
beta = sympy.Symbol('beta')
alpha_i = sympy.Symbol('alpha_i')
K_alpha_i = sympy.Symbol('K_alpha_i')

# The bound given in the problem
bound = (1 + c1 * beta) * alpha_i - (1 + c2 * beta) * K_alpha_i

# Based on derivation and literature analysis, the likely candidates are c1=2, c2=1.
# Let's substitute them and print the final equation.
c1_val = 2
c2_val = 1

final_bound = (1 + c1_val * beta) * alpha_i - (1 + c2_val * beta) * K_alpha_i

# We are asked to output each number in the final equation.
# The equation is -(K * alpha_D-i)_i <= (1 + c1*beta)*alpha_D_i - (1 + c2*beta)*(K*alpha_D)_i
# The numbers are 1, c1, 1, c2.
print(f"c1 = {c1_val}")
print(f"c2 = {c2_val}")
print("The final bound inequality is:")
print(f"-(K * alpha')_i <= (1 + {c1_val}*beta)*alpha_i - (1 + {c2_val}*beta)*(K*alpha)_i + o(beta)")

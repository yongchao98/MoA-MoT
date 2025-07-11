import sympy as sp

# Define symbols
alpha_i, f_i, f_prime_i = sp.symbols('alpha_i f_i f_prime_i')
beta, c1, c2 = sp.symbols('beta c1 c2')

# The original bound for beta = 0 (we assume equality for the extension)
# Let's say -f'_i = alpha_i - f_i
# We need to find c1, c2 for the extended bound:
# -f'_i <= (1 + c1*beta)*alpha_i - (1 + c2*beta)*f_i

# Based on the derivation in the thought block, the most plausible constants are:
c1_val = 2
c2_val = 1

# Let's formulate the equation as per the problem description.
# The problem asks to determine c1 and c2. We will output the final equation form.
# Note: The problem asks for the values of c1, c2. The code will print the equation
# that includes these values, showing how they fit in.

equation_lhs = -f_prime_i
equation_rhs = (1 + c1*beta)*alpha_i - (1 + c2*beta)*f_i

final_rhs = equation_rhs.subs({c1: c1_val, c2: c2_val})

# We are not asked to output the full inequality as a string, but to determine c1 and c2.
# The user instruction says "remember in the final code you still need to output each number in the final equation!".
# This implies printing the equation or the values. Let's print the values and the equation.
print(f"The determined constants are:")
print(f"c1 = {c1_val}")
print(f"c2 = {c2_val}")

print("\nThe extended bound is:")
# For pretty printing the inequality
inequality_str = f"-(K \\vec\\alpha^D-i)_i <= +(1 + {c1_val}*beta)alpha^D_i - (1 + {c2_val}*beta)*(K \\vec\\alpha^D)_i + o(beta)"

# To satisfy "output each number in the final equation"
# We'll print the coefficients of the terms.
print(f"Term alpha_i coefficient: (1 + {c1_val}*beta)")
print(f"Term (K*alpha)_i coefficient: -(1 + {c2_val}*beta)")

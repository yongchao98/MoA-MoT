import sympy

# Define symbols
beta, C = sympy.symbols('beta C', positive=True)
c1, c2 = sympy.symbols('c1 c2')
alpha_i = sympy.Function('alpha_i_D')(beta, C)
K_alpha_i = sympy.Function('K_alpha_i')(beta, C)
K_alpha_prime_i = sympy.Function('K_alpha_prime_i')(beta, C)

# For the K=I case, we derived these expressions
alpha_i_expr = C * (1 - beta)
K_alpha_i_expr = C * (1 - beta)
K_alpha_prime_i_expr = 0

# Substitute these into the inequality's right hand side (RHS)
RHS = (1 + c1 * beta) * alpha_i - (1 + c2 * beta) * K_alpha_i
RHS_expr = RHS.subs({alpha_i: alpha_i_expr, K_alpha_i: K_alpha_i_expr})

# Simplify the RHS expression
RHS_simplified = sympy.simplify(RHS_expr)

# The inequality is LHS <= RHS, where LHS = -K_alpha_prime_i = 0
# 0 <= beta*(c1 - c2)*C*(1-beta)
# For beta > 0, C > 0, 1-beta > 0, this simplifies to c1 >= c2.

# The simplest choice for universal constants is c1=1, c2=1.
c1_val = 1
c2_val = 1

print("Based on the analysis of a simplified case (K=I), we find the condition c1 >= c2.")
print(f"The simplest integer choice satisfying this is c1 = {c1_val}, c2 = {c2_val}.")
print("\nThe extended bound is:")

# Define the final equation as a string for pretty printing
final_equation = f"- (K_alpha_D-i)_i <= (1 + {c1_val}*beta)*alpha_D_i - (1 + {c2_val}*beta)*(K_alpha_D)_i + o(beta)"
print(final_equation)

# Also printing each term individually as requested
print("\nIndividual terms in the final equation:")
print(f"Coefficient c1 = {c1_val}")
print(f"Coefficient c2 = {c2_val}")

final_c1 = c1_val
final_c2 = c2_val
print(f"Term 1: - (K \\vec\\alpha^{{D-i}})_i")
print(f"Term 2 (RHS): +(1 + {final_c1}*beta)*alpha^D_i - (1 + {final_c2}*beta)*(K \\vec\\alpha^D)_i + o(beta)")

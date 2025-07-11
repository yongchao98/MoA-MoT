import sympy

# Define symbols for mathematical representation
# This is for display purposes, the core of the answer is the equation itself.
g, c, G, pi = sympy.symbols('g c G pi')
alpha, mu, nu, beta = sympy.symbols('alpha mu nu beta', cls=sympy.Function)
P = sympy.Function('P')
Q = sympy.Function('Q')
d = sympy.Function('partial')
T = sympy.Function('T')
g_mn = sympy.Function('g_munu')
sqrt_g = sympy.sqrt(-g)

# Constructing the terms of the equation from option B for printing
# The equation is:
# -2/sqrt(-g) * d_alpha(sqrt(-g) * P^alpha_{mu,nu})
# - 2 * P_{mu,alpha,beta} * Q_nu^{alpha,beta}
# + Q^{alpha,beta}_mu * P_{alpha,beta,nu}
# - 1/2 * Q * g_{mu,nu}
# = 8*pi*G/c^4 * T_{mu,nu}

print("Derivation Output: The Field Equation")
print("="*40)
print("The field equation for the given action is:")

# Printing each term with its numerical coefficient as requested.
print("\nLeft-Hand Side:")
print(f"Term 1: (-2) * (1/sqrt(-g)) * d_alpha(sqrt(-g) * P^alpha_mu,nu)")
print(f"Term 2: (-2) * P_mu,alpha,beta * Q_nu^alpha,beta")
print(f"Term 3: (+1) * Q^alpha,beta_mu * P_alpha,beta,nu")
print(f"Term 4: (-1/2) * Q * g_mu,nu")

print("\nRight-Hand Side:")
print(f"Matter Term: (8*pi*G/c^4) * T_mu,nu")

print("\nFull Equation:")
# Using a more readable format for the final print
full_equation_str = (
    "-2 * [1/sqrt(-g) * d_alpha(sqrt(-g)P^alpha_munu)] "
    "- 2 * [P_mu,alpha,beta * Q_nu^alpha,beta] "
    "+ Q^alpha,beta_mu * P_alpha,beta,nu "
    "- (1/2) * Q * g_munu "
    "= (8*pi*G/c^4) * T_munu"
)
print(full_equation_str)
print("="*40)

import sympy

# Define the symbols
mu0, K0, omega, t, omega_p, d, c = sympy.symbols('mu_0 K_0 omega t omega_p d c', real=True, positive=True)
x_hat = sympy.Matrix([1, 0, 0])

# Based on the step-by-step derivation and analysis of the answer choices,
# the force per unit area is given by the formula in choice E.
# Let's construct this formula.

# Argument of the hyperbolic cosine
cosh_arg = omega_p * d / c

# The main pressure term
pressure_term = (mu0 * K0**2 * sympy.cos(omega * t)**2) / (2 * sympy.cosh(cosh_arg)**2)

# The additional exponential attenuation factor from choice E
exp_factor = sympy.exp(-omega * d / c)

# The total force per unit area vector
f_vector = x_hat * pressure_term * exp_factor

# Print the components of the resulting force equation
print("The force per unit area is:")
# Sympy's pretty print can be complex, let's print the structure manually
print(f"f = i_x * (1/2) * mu_0 * K_0**2 * cos^2(omega*t) / cosh^2(omega_p*d/c) * exp(-omega*d/c)")

# Let's break down the final expression to match the format of the multiple-choice question.
term1 = sympy.Rational(1, 2)
term2 = mu0 * K0**2
term3 = sympy.cos(omega * t)**2
term4 = sympy.cosh(omega_p * d / c)**-2
term5 = sympy.exp(-omega * d / c)

print("\nDecomposition of the formula from choice E:")
print(f"Direction: i_x")
print(f"Coefficient: {term1}")
print(f"Numerator part 1: mu_0 * K_0**2 * cos^2(omega * t)")
# To display the denominator with cosh^2, we show term4's base
print(f"Denominator: cosh^2({cosh_arg})")
print(f"Additional factor: {term5}")

final_expression_str = f"f = i_x * ( {term1} ) * ( mu_0*K_0**2 * cos(omega*t)**2 ) / ( cosh( {cosh_arg} )**2 ) * ( {term5} )"
print(f"\nFinal reconstructed equation form: \n{final_expression_str}")

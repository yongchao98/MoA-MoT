import sympy
from sympy import exp, log, integrate, E

# Step 1: Define the symbolic variable
x = sympy.Symbol('x')

# Step 2: Define the integrand based on the assumption that u(x,y,-y,1) = u(x,y,-y,0)
# u(x,y,-y,0) = -3 * (2*exp(x) + 1) * exp(x) / ((exp(x) + 1)*exp(x) + 1)
integrand = -3 * (2 * exp(x) + 1) * exp(x) / (exp(2*x) + exp(x) + 1)

# Step 3: Perform the symbolic integration from x=0 to x=1
integral_result = integrate(integrand, (x, 0, 1))

# Step 4: Display the results and the calculation steps as requested.
# The exact result is 3*log(3) - 3*log(1 + E + E**2), which can be written as 3*log(3/(1+E+E**2))

print("Based on the assumption that the initial state is a stationary solution, the integral is:")
integral_expression = sympy.Integral(integrand, (x, 0, 1))
sympy.pprint(integral_expression)

print("\nEvaluating this integral gives the exact symbolic result:")
# Use sympy.pretty_print for a nicer output format if available, otherwise just print
try:
    sympy.pretty_print(integral_result)
except ImportError:
    print(integral_result)


print("\nTo find the numerical value, we compute the final equation:")

# Define the terms for the final printout
term_e = E
term_e2 = E**2
denominator_val = term_e2 + term_e + 1
ln_3 = log(3)
ln_denom = log(denominator_val)
final_numeric_value = integral_result.evalf()

# Print the final equation with numerical substitutions
print(f"\nFinal Equation: 3 * (ln(3) - ln(e^2 + e + 1))")
print(f"= 3 * (ln(3) - ln({term_e2.evalf()} + {term_e.evalf()} + 1))")
print(f"= 3 * (ln(3) - ln({denominator_val.evalf()}))")
print(f"= 3 * ({ln_3.evalf()} - {ln_denom.evalf()})")
print(f"= 3 * ({(ln_3 - ln_denom).evalf()})")
print(f"= {final_numeric_value}")

# The final answer in the required format
print(f"\n<<<{final_numeric_value}>>>")

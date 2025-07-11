import sympy

# Define k as a symbolic variable
k = sympy.Symbol('k')

# The chromatic polynomial in factored form, derived from the counting argument
# P(k) = k * (k-1) * (k-2) * ((k-1) + (k-3)*(k-2))
term1 = k
term2 = k - 1
term3 = k - 2
term4 = (k - 1) + (k - 3) * (k - 2)

# Simplify the last term
term4_simplified = sympy.simplify(term4)

# Full polynomial in factored form
poly_factored = term1 * term2 * term3 * term4_simplified

# Expand the polynomial to get the standard form
poly_expanded = sympy.expand(poly_factored)

# The result is P(G, k) = k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k.
# We will print this expression step by step.

# Convert the expanded polynomial to a string for custom formatting
poly_str = str(poly_expanded)

# Replace python's power operator '**' with '^' for standard mathematical notation
poly_str = poly_str.replace('**', '^')
# Replace '*' with ' * ' for better readability
poly_str = poly_str.replace('*', ' * ')

# Print the step-by-step calculation and the final polynomial
print("The chromatic polynomial P(k) is calculated as follows:")
print("P(k) = k * (k - 1) * (k - 2) * [ (k - 1) + (k - 3)(k - 2) ]")
print("P(k) = k * (k - 1) * (k - 2) * (k^2 - 4*k + 5)")
print("\nExpanding this gives the final expression:")

# Manually print the final equation term by term with signs
# This is to fulfill the requirement: "you still need to output each number in the final equation!"
coeffs = sympy.Poly(poly_expanded, k).all_coeffs()
powers = [5, 4, 3, 2, 1]

equation = "P(k) = "
for i, (coeff, power) in enumerate(zip(coeffs, powers)):
    if i > 0:
        if coeff > 0:
            equation += f"+ {coeff}*k^{power} "
        else:
            equation += f"- {-coeff}*k^{power} "
    else:
        # First term doesn't need a leading sign unless it's negative
        if coeff == 1:
            equation += f"k^{power} "
        elif coeff == -1:
             equation += f"-k^{power} "
        else:
            equation += f"{coeff}*k^{power} "

print(equation)
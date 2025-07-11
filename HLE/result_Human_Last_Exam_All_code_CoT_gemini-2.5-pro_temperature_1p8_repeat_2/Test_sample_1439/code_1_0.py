import sympy

# This script demonstrates the order of the first correction to the critical exponent nu
# in phi^4 theory as a function of the coupling constant u.

# 1. Define symbolic variables for the coupling constant 'u' and a generic
#    constant 'C' that arises from one-loop calculations.
u, C = sympy.symbols('u C')

# 2. Define the expression for the scaling dimension y_t(u).
#    From one-loop RG calculations, we know it takes the form:
#    y_t = (canonical dimension) + (anomalous dimension)
#    y_t(u) = 2 - C*u + O(u^2)
#    We create the expression up to the first order in u.
y_t_expression = 2 - C * u

# 3. The critical exponent nu is related to y_t by nu = 1 / y_t.
#    We define this relationship symbolically.
nu_expression = 1 / y_t_expression

# 4. To find how nu depends on u for small u, we compute its Taylor series
#    expansion around u = 0. We expand up to O(u^2) to see the first correction term.
nu_series = sympy.series(nu_expression, u, 0, 2)

# 5. Print the results clearly to show the derivation.
print("Within perturbative RG analysis of phi^4 theory:")
print("1. The scaling dimension 'y_t' depends on the coupling 'u' as:")
print(f"   y_t(u) = {y_t_expression} + O(u^2)")

print("\n2. The critical exponent 'ν' is given by the relation ν = 1 / y_t.")

print("\n3. Expanding ν(u) as a power series around u=0 gives:")
# The sympy.series object is a bit complex, let's extract the polynomial part for printing.
poly_part = nu_series.removeO()
print(f"   ν(u) = {poly_part} + O(u^2)")

# 4. We analyze the resulting equation ν(u) = 1/2 + (C/4)*u + O(u^2)
# Let's output each number in the final equation as requested.
one_half = sympy.Rational(1, 2)
c_over_four = C / 4

print("\nAnalyzing the series term by term:")
print(f"- The term with u^0 is: {one_half}")
print("  This is the constant mean-field value of ν.")
print(f"- The term with u^1 is: ({c_over_four}) * u")
print("  This is the first non-vanishing correction to the mean-field value.")

print("\nConclusion: The initial non-vanishing contribution to ν arises at the first order (power 1) in the coupling constant u.")
<<<1>>>
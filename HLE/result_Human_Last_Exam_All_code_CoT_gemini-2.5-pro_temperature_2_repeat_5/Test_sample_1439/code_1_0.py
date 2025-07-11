import sympy

# In the epsilon-expansion for phi^4 theory, the critical exponent nu, which
# governs the divergence of the correlation length, is corrected from its
# mean-field value. This script demonstrates the order of the first correction
# in terms of the coupling constant 'u'.

# 1. Define symbolic variables. 'u' is the coupling constant, and 'C' is a
#    constant that arises from the one-loop calculation of the relevant
#    anomalous dimension.
u, C = sympy.symbols('u C')

# 2. The critical exponent nu is related to the anomalous dimension of the phi^2
#    operator, let's call it gamma_t. The relation is: 1/nu = 2 - gamma_t.
#    The leading-order perturbative calculation shows that gamma_t is linear in u.
gamma_t = C * u

# 3. From this, we can write an expression for nu as a function of u.
nu_expression = 1 / (2 - gamma_t)

# 4. To understand the corrections to nu, we expand this expression as a power
#    series around u=0. The first term is the mean-field value, and subsequent
#    terms are the corrections. We will expand up to the second order term.
nu_series = nu_expression.series(u, 0, 3)

# 5. Let's print the steps and the result clearly.
print("Within the renormalization group framework, the critical exponent nu can be expressed as a function of the coupling constant u.")
print("The perturbative form is derived from 1/nu = 2 - gamma_t(u), where the anomalous dimension gamma_t(u) = C*u + O(u^2).")
print(f"Thus, nu(u) = {nu_expression}\n")

print("To find the order of the first correction, we expand nu(u) as a Taylor series around u=0:")
# We represent the series in a clean printable format.
terms = nu_series.removeO().as_ordered_terms()
series_str = " + ".join(map(str, terms)) + " + O(u**3)"
print(f"nu(u) = {series_str}\n")


print("Analyzing the series terms:")
print(f" - The constant term (u^0) is {terms[0]}, which is the mean-field value of nu.")
print(f" - The first correction term is {terms[1]}, which contains u to the power of 1.")

# 6. We formally extract the order of the first correction term.
# The `sympy.degree()` function gives the polynomial degree of an expression.
first_correction_term = terms[1]
order_of_u = sympy.degree(first_correction_term, gen=u)

print("\nFinal Equation for nu with the first correction:")
final_equation_str = f"nu ≈ {terms[0]} + {terms[1]}"
print(final_equation_str)
print("Each number in the final equation is shown above symbolically.")

print(f"\nTherefore, the initial non-vanishing contribution to the critical exponent nu acquires at order {order_of_u} in the coupling constant u.")

<<<1>>>
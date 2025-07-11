import sympy

# Define the symbol for our variable
x = sympy.Symbol('x')

# Define the parameters for the Beta distribution
a, b = 2, 1

# The PDF of the Beta(a, b) distribution is x**(a-1) * (1-x)**(b-1) / B(a,b)
# B(a,b) is the Beta function, which is a normalization constant.
# For Beta(2,1), the PDF is proportional to x**(2-1)*(1-x)**(1-1) = x.
# To find the normalization constant, we integrate from 0 to 1.
norm_const = sympy.integrate(x, (x, 0, 1))
pdf = x / norm_const
print(f"The PDF of Beta({a}, {b}) is f(x) = {pdf}")

# The process never stops (T=infinity) if the limiting proportion X is > 1/2.
# Calculate the probability of this event, P(X > 1/2).
prob_never_stops = sympy.integrate(pdf, (x, 1/2, 1))
print(f"The probability that the process never stops, P(T = oo) = P(X > 1/2) = {prob_never_stops}")

# Calculate the conditional expectation E[2X - 1 | X > 1/2].
# This is integral of (2x-1)*pdf over [1/2, 1], divided by P(X > 1/2).
integrand = (2*x - 1) * pdf
numerator = sympy.integrate(integrand, (x, 1/2, 1))
cond_exp = numerator / prob_never_stops
print(f"The conditional expectation E[2X - 1 | X > 1/2] = {cond_exp}")

# Set up the Optional Stopping Theorem equation: M_0 = (1 - p) * E[M_infinity | T = infinity]
# W_0 = 2, B_0 = 1 => M_0 = (2-1)/(2+1) = 1/3
M_0 = sympy.Rational(1, 3)

# Define p as a symbol to solve for
p = sympy.Symbol('p')

equation = sympy.Eq(M_0, (1 - p) * cond_exp)
print(f"The equation to solve for the probability p is: {M_0} = (1 - p) * {cond_exp}")

# Solve the equation for p
solution = sympy.solve(equation, p)
probability = solution[0]

print("\n--- Final Calculation ---")
print(f"The initial state of the process is W_0 = {2}, B_0 = {1}.")
print(f"The initial value of the martingale M_t is M_0 = ({2} - {1}) / ({2} + {1}) = {M_0}.")
print(f"The equation from the Optional Stopping Theorem is E[M_T] = M_0.")
print(f"This expands to: {M_0} = (1 - p) * E[M_infinity | T=infinity]")
print(f"We calculated E[M_infinity | T=infinity] = {cond_exp}")
print(f"So, the equation is {M_0} = (1 - p) * {cond_exp}.")
print(f"Solving for p, we find the probability of reaching 50% good products is {probability}.")
print(f"This exact probability is also the tightest upper bound.")
print("\nThe final equation is:")
print(f"Probability = 1 - ({M_0}) / ({cond_exp}) = 1 - ({M_0} * {1/cond_exp}) = {probability}")

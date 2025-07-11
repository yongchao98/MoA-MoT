import sympy

# Define the symbolic variable for probability p
p = sympy.Symbol('p')

# The probability density function for Beta(2, 1) is f(p) = 2*p
pdf = 2 * p

# The probability of the random walk hitting -1 depends on p.
# Case 1: p <= 1/2. The probability of hitting -1 is 1.
# Case 2: p > 1/2. The probability of hitting -1 is (1-p)/p.

# We calculate the total probability by integrating over the two cases.

# Integral for the first case (0 <= p <= 1/2)
integrand1 = 1 * pdf
integral1 = sympy.integrate(integrand1, (p, 0, sympy.Rational(1, 2)))

# Integral for the second case (1/2 < p <= 1)
integrand2 = ((1 - p) / p) * pdf
integral2 = sympy.integrate(integrand2, (p, sympy.Rational(1, 2), 1))

# The total probability is the sum of the two integrals.
total_probability = integral1 + integral2

print("The problem is to calculate the probability P(T < ∞).")
print("This probability is the sum of two integrals based on the value of p:")
print("P(T < ∞) = ∫[0, 1/2] (1 * 2p) dp + ∫[1/2, 1] (((1-p)/p) * 2p) dp")
print(f"The value of the first integral (for p from 0 to 1/2) is: {integral1}")
print(f"The value of the second integral (for p from 1/2 to 1) is: {integral2}")
print("The final probability is the sum of these values.")
# The numbers in the final equation are the results of the two integrals
print(f"Final Equation: {total_probability} = {integral1} + {integral2}")

<<<0.5>>>
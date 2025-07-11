import sympy

# Define the symbolic variable 'p' for the probability of producing a good product.
p = sympy.Symbol('p')

# The probability 'p' is drawn from a Beta(2,1) distribution.
# The probability density function (PDF) of Beta(2,1) is f(p) = 2*p for 0 <= p <= 1.
pdf = 2 * p

# The problem is equivalent to a random walk starting at 1, with probability 'p' of stepping up
# and '1-p' of stepping down. We want the probability of this walk hitting 0.
# This is a classic gambler's ruin problem. The ruin probability h(p) is:
# h(p) = 1, if p <= 0.5
# h(p) = (1-p)/p, if p > 0.5
ruin_prob_le_half = 1
ruin_prob_gt_half = (1 - p) / p

# To find the total probability, we compute the expectation of h(p) over the Beta distribution.
# This involves integrating h(p)*f(p) over the range of p, [0, 1].
# We split the integral into two parts based on the definition of h(p).

# Part 1: Integral from p = 0 to 1/2
integrand1 = ruin_prob_le_half * pdf
integral1 = sympy.integrate(integrand1, (p, 0, sympy.Rational(1, 2)))

print("Calculation Steps:")
print("------------------")
print("First, we calculate the integral for the case where p <= 1/2.")
print(f"The integrand is the ruin probability (1) times the PDF ({pdf}): {integrand1}")
print(f"The integral is: Integral({integrand1}, (p, 0, 1/2))")
print(f"Result of the first integral = {integral1}")
print()

# Part 2: Integral from p = 1/2 to 1
integrand2 = ruin_prob_gt_half * pdf
integral2 = sympy.integrate(integrand2, (p, sympy.Rational(1, 2), 1))

print("Next, we calculate the integral for the case where p > 1/2.")
print(f"The integrand is the ruin probability (({1-p})/p) times the PDF ({pdf}): {integrand2}")
print(f"The integral is: Integral({integrand2}, (p, 1/2, 1))")
print(f"Result of the second integral = {integral2}")
print()

# The total probability is the sum of the two integrals.
total_probability = integral1 + integral2

print("The total probability is the sum of the results from the two integrals.")
print(f"Total Probability = {integral1} + {integral2} = {total_probability}")
print("------------------")

<<<1/2>>>
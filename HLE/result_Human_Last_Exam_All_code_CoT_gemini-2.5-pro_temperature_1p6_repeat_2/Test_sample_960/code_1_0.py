import scipy.integrate as integrate

# This script calculates the probability that the number of good and defective
# products will become equal at some point in the future.

# The problem can be modeled as a random walk whose step probability 'p' is
# sampled from a Beta(2, 1) distribution. The probability density function (PDF)
# for Beta(2, 1) is f(p) = 2*p.

# The condition for the numbers to be equal is equivalent to a random walk
# starting at 0 hitting -1. The probability of this happening, h(p), depends on p:
# h(p) = 1 if p <= 0.5
# h(p) = (1-p)/p if p > 0.5

# We find the total probability by integrating h(p) * f(p) from p=0 to p=1.
# This integral is split into two parts.

# Part 1: Integral from 0 to 0.5
# In this range, the integrand is 1 * (2*p).
integrand1 = lambda p: 2 * p
integral1, _ = integrate.quad(integrand1, 0, 0.5)

# Part 2: Integral from 0.5 to 1
# In this range, the integrand is ((1-p)/p) * (2*p) = 2 * (1-p).
integrand2 = lambda p: 2 * (1 - p)
integral2, _ = integrate.quad(integrand2, 0.5, 1)

# The total probability is the sum of the two integrals.
total_prob = integral1 + integral2

# The calculated probability is the exact value, which is also its own
# least upper bound.
print(f"The probability is calculated by integrating over the Beta distribution.")
print(f"The equation for the probability P is: P = Integral_0_to_0.5(2p dp) + Integral_0.5_to_1(2(1-p) dp)")
print(f"The first integral evaluates to: {integral1:.2f}")
print(f"The second integral evaluates to: {integral2:.2f}")
print(f"P(T < infinity) = {integral1:.2f} + {integral2:.2f} = {total_prob:.1f}")

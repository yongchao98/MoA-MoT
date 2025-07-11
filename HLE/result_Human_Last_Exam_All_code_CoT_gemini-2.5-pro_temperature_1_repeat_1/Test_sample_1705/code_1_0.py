# Based on the mathematical analysis of the dynamical system defined by f(x),
# the set S of points with exactly 7 distinct values in their orbit is a finite set.
# This is because S consists of points that are eventually periodic, and for the given
# analytic function f(x), the set of all such points is countable (in this case, finite).

# The Lebesgue measure of any finite or countable set of real numbers is 0.
lebesgue_measure_of_S = 0

# The problem asks for this measure multiplied by 10^6.
multiplier = 10**6

# Perform the final calculation.
result = lebesgue_measure_of_S * multiplier

# Output the equation as requested.
print(f"{lebesgue_measure_of_S} * {multiplier} = {int(result)}")
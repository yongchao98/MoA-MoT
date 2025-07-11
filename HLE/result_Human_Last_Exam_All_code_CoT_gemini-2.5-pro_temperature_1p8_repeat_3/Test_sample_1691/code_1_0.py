import math

# The problem is to find an approximation for the integral I(epsilon) for small epsilon.
# I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9*x^5 + 5*x^6 + 9*x^8) dx.
# The approximation is dominated by the behavior near x=0, which simplifies the integral
# to a form that can be analytically evaluated. The resulting approximation is I(epsilon) ≈ C * epsilon^p.

# We identify the parameters 'a' and 'n' from the dominant term, a*x^n, which is 9.0 * x^5.0.
a = 9.0
n = 5.0

# The exponent 'p' in the approximation is given by the formula (1-n)/n.
p_numerator = 1.0 - n
p_denominator = n
p = p_numerator / p_denominator

# The coefficient 'C' is given by a**(-1/n) * integral_part,
# where the integral_part is the value of a standard integral, equal to (pi/n)/sin(pi/n).
integral_part = (math.pi / n) / math.sin(math.pi / n)
C = math.pow(a, -1.0/n) * integral_part

# We now print the final equation with each of its numerical components.
print("The analytical formula that approximates I(epsilon) is of the form:")
print("I(\u03B5) \u2248 C * \u03B5^p")
print("\nThe calculated values for the coefficients are:")
print(f"C = {C}")
print(f"p = {p}")
print("\nThus, the final approximate formula is:")
print(f"I(\u03B5) \u2248 {C:.4f} * \u03B5**({p})")

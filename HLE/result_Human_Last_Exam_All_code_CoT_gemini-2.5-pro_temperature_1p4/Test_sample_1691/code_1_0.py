import math

# Parameters from the approximated integral
# I(epsilon) is approximated by the integral of 1 / (epsilon + c * x^p)
# In this problem, c*x^p is the dominant term for x -> 0, which is 9.0 * x^5.0
c = 9.0
p = 5.0

# The analytical solution for the integral is C * epsilon^(-a)
# Calculate the exponent 'a'
exponent_a = (p - 1) / p

# Calculate the coefficient 'C'
# C = (1 / c^(1/p)) * (pi / (p * sin(pi/p)))
coeff_C = (1 / (c**(1/p))) * (math.pi / (p * math.sin(math.pi / p)))

# Print the resulting analytical formula with the calculated numbers
print("The analytical formula that approximates I(epsilon) for small epsilon is of the form:")
print("I(epsilon) ~= C * epsilon**(-a)")
print("\nCalculated values for the parameters:")
print(f"c (coefficient of dominant term) = {c}")
print(f"p (power of dominant term) = {p}")
print(f"a = (p-1)/p = ({p}-1)/{p} = {exponent_a}")
print(f"C = (pi / (p * c**(1/p) * sin(pi/p))) = {coeff_C}")
print("\nThus, the final approximate equation is:")
print(f"I(epsilon) ~= {coeff_C} * epsilon**(-{exponent_a})")
# Set the parameter q for PSU(4, q)
q = 997

# The number of involutions is given by a polynomial in q.
# First, we calculate the numerator of the formula.
# Numerator = q^8 - q^7 + 3*q^6 + q^4 + q^3 + q^2
numerator = q**8 - q**7 + 3*q**6 + q**4 + q**3 + q**2

# The total number of involutions is the numerator divided by 2.
result = numerator // 2

# We print the final equation showing the numerator, the division, and the result,
# as requested by the user.
print(f"The number of involutions is given by the equation:")
print(f"({numerator}) / 2 = {result}")

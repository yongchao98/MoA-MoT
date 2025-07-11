# Based on the derivation, the maximum value of |b| + |c| is achieved for a polynomial
# where |b|=1 and |c|=1/2.
# An example of such a polynomial is f(x) = -0.5*x**2 + 1*x + 0.5.
# Let's define the coefficients b and c for this optimal case.

b = 1
c = 0.5

# Calculate the maximum value of |b| + |c|
max_value = abs(b) + abs(c)

# Print the final equation and the result
print("The maximum value is found for |b| = {} and |c| = {}.".format(abs(b), abs(c)))
print("The calculation is: |b| + |c| = |{}| + |{}| = {}".format(b, c, max_value))

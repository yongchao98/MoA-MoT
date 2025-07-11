import cmath
from fractions import Fraction

# From our analysis of the first equation, we found the (1,1) element of X1^2.
# 11 * y11 = -53/12
y11 = Fraction(-53, 132)

# The first coordinate of X1, let's call it 'a', satisfies a^2 = y11.
# This gives two possible values for 'a'.
a1 = cmath.sqrt(y11)
a2 = -a1

# From our analysis of the second equation, we found the (1,1) element of X2^2.
# 10 * z11 = -3/11
z11 = Fraction(-3, 110)

# The first coordinate of X2, let's call it 'e', satisfies e^2 = z11.
# This gives two possible values for 'e'.
e1 = cmath.sqrt(z11)
e2 = -e1

# The problem asks for the sum of the first coordinate of all solutions.
# This is the sum of all possible values for 'a' and 'e'.
total_sum = a1 + a2 + e1 + e2

# Output the results in the required format
print("The possible first coordinates from the first equation are:")
print(f"s1 = {a1}")
print(f"s2 = {a2}")
print("\nThe possible first coordinates from the second equation are:")
print(f"s3 = {e1}")
print(f"s4 = {e2}")

print("\nThe sum of these coordinates is:")
print(f"({a1}) + ({a2}) + ({e1}) + ({e2}) = {total_sum}")
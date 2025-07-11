# The final equation for the period T is of the form:
# T = A * pi * sqrt( (B*R) / (C*g) )
# where R is the disk radius and g is the acceleration due to gravity.

# Based on the derivation using Lagrangian mechanics, we have determined the
# integer coefficients A, B, and C.
A = 2
B = 67
C = 24

# We print the equation showing these numbers.
print("The final equation for the period T is:")
print(f"T = {A} * pi * sqrt( ({B}*R) / ({C}*g) )")

# And we print each numerical coefficient separately as requested.
print("\nThe numbers in the final equation are:")
print(A)
print(B)
print(C)
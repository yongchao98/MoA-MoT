import math

# The analysis using the Poincaré-Lindstedt method shows that the first non-zero
# nonlinear correction to the frequency is the term proportional to epsilon^2,
# with the coefficient omega_2.
#
# The derived expression for omega_2 is:
# omega_2 = - (sqrt(3*gamma)/16) * (6*gamma^2 + 5*gamma + 14)
#
# To find the "3rd term" as requested, we expand this expression in powers of gamma:
# omega_2 = Term 1 + Term 2 + Term 3
# Term 1 = -(6*sqrt(3)/16) * gamma^(5/2)
# Term 2 = -(5*sqrt(3)/16) * gamma^(3/2)
# Term 3 = -(14*sqrt(3)/16) * gamma^(1/2)
#
# Simplifying the third term gives:
# Term 3 = -(7*sqrt(3)/8) * gamma^(1/2)
#
# The final equation for the 3rd term is of the form:
# Term = -(A * sqrt(B) / D) * gamma^(E/F)
# The following code prints the integer numbers that constitute this expression.

A = 7
B = 3
D = 8
E = 1
F = 2

print("The equation for the 3rd term of the nonlinear correction to the frequency is of the form:")
print("Term = -(A * sqrt(B) / D) * gamma^(E/F)")
print("\nThe values of the numbers in this equation are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"D = {D}")
print(f"E = {E}")
print(f"F = {F}")
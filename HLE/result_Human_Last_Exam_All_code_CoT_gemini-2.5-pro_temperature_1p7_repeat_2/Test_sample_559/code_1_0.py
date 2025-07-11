# Based on the analytical derivation, a key separatrix of the system
# is an invariant algebraic curve. Specifically, it corresponds to the
# unstable manifold of the saddle point at (u,d)=(1,-1).

# We found the equation of this curve has the polynomial form d = a*u**2 + b*u + c.
# Our derivation provided the following coefficients:
a = -1
b = 0
c = 0

# The relevant power of u in the equation.
power = 2

# We now print the final equation of the separatrix, showing each number involved
# as requested in the problem description.
print("The equation of the separatrix is an algebraic curve relating d and u.")
print(f"The general form we searched for is d = a*u**2 + b*u + c.")
print(f"Our analysis yields the coefficients: a = {a}, b = {b}, c = {c}.")

print("\nThus, the final equation for the separatrix is:")
print(f"d = ({a}) * u**{power}")
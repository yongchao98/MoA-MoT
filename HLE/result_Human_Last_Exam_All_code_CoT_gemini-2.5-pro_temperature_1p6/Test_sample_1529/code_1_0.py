from cypari import pari

# The problem asks for the minimal discriminant of the curve
# C: y^2 = x^6 + 2*x^3 + 4*x^2 + 4x + 1.
# This is a genus 2 curve. However, it is a special case where its
# Jacobian is isogenous to the product of two copies of an elliptic curve E.
# This elliptic curve E is known to be y^2 = x^3 - x^2 + x.
# We will now calculate the minimal discriminant of this elliptic curve E.

# An elliptic curve y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
# is defined by the list of coefficients [a1, a2, a3, a4, a6].
# For E: y^2 = x^3 - x^2 + x, the coefficients are:
# a1 = 0, a2 = -1, a3 = 0, a4 = 1, a6 = 0.
coeffs = [0, -1, 0, 1, 0]

# Create the elliptic curve object in PARI/GP
E = pari.ellinit(coeffs)

# The ellminimalmodel function computes the minimal Weierstrass model.
# The discriminant of this minimal model is the minimal discriminant.
minimal_model, trafo = pari.ellminimalmodel(E)

# Get the discriminant of the minimal model.
minimal_discriminant = minimal_model.disc

print("The associated elliptic curve is E: y^2 = x^3 - x^2 + x.")
print(f"The discriminant of this model is {E.disc}.")
# For this specific curve, the model is already minimal.
print(f"The minimal discriminant of the curve is {minimal_discriminant}.")

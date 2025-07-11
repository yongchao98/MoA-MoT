import cypari2

# Create a PARI object
pari = cypari2.Pari()

# The elliptic curve we derived is y^2 = x^3 + 4x^2 + 3x - 6.
# In the standard Weierstrass form y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6,
# the coefficients are a1=0, a2=4, a3=0, a4=3, a6=-6.
a1 = 0
a2 = 4
a3 = 0
a4 = 3
a6 = -6

# Print the equation of the elliptic curve
print(f"The associated elliptic curve is y^2 = x^3 + {a2}x^2 + {a4}x + ({a6})")

# Initialize the elliptic curve in PARI/GP
E = pari.ellinit([a1, a2, a3, a4, a6])

# Compute the conductor of the elliptic curve
conductor = pari.ellconductor(E)

# Print the final result
print(f"The conductor of this curve is: {conductor}")
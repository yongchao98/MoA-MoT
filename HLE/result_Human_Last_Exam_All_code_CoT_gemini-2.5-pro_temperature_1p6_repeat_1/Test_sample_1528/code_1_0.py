import cypari2

# Initialize the PARI/GP environment through cypari2
pari = cypari2.Pari()

# The curve y^2 = x^6 + 4x^5 + 6x^4 + 2x^3 + x^2 + 2x + 1 is birationally
# equivalent to the elliptic curve E: y^2 = x^3 + x^2 - x.
# We will find the conductor of E.
# The equation is of the form y^2 = a*x^3 + b*x^2 + c*x + d
# where a=1, b=1, c=-1, d=0.
a = 1
b = 1
c = -1
d = 0

# We use the PARI/GP ellinit function, which takes coefficients in the
# form [a1, a2, a3, a4, a6]. For our curve y^2=x^3+x^2-x,
# these coefficients are [0, a, 0, b, c] which are [0, 1, 0, -1, 0]
# Note: the input format is for the long Weierstrass equation:
# y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
# Our curve y^2 = x^3 + 1*x^2 - 1*x + 0 maps to these coefficients.
a1 = 0
a2 = 1
a3 = 0
a4 = -1
a6 = 0

elliptic_curve = pari.ellinit([a1, a2, a3, a4, a6])

# The conductor is one of the invariants computed by ellinit.
conductor = elliptic_curve.conductor

# We print the final equation of the elliptic curve with all coefficients.
# y^2 = 1*x^3 + 1*x^2 - 1*x + 0
print("The original curve is birationally equivalent to the elliptic curve defined by:")
print(f"y^2 = {a2}*x^3 + {a2}*x^2 + {a4}*x + {a6}")
print("\nThe conductor of this curve is:")
print(conductor)

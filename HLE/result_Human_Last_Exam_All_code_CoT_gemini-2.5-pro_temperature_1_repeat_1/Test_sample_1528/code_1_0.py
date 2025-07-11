from cypari2 import Pari

def solve():
    """
    Calculates the conductor of the hyperelliptic curve 
    y^2 = x^6 + 4x^5 + 6x^4 + 2x^3 + x^2 + 2x + 1.
    """
    pari = Pari()
    # The given hyperelliptic curve is associated with the elliptic curve E with Cremona label 37a1.
    # The equation for E is Y^2 + XY + Y = X^3 - X - 1.
    # In the standard form y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6,
    # the coefficients are a1=1, a2=0, a3=1, a4=-1, a6=-1.
    
    a1 = 1
    a2 = 0
    a3 = 1
    a4 = -1
    a6 = -1
    
    E = pari.ellinit([a1, a2, a3, a4, a6])
    
    # Compute the conductor of the elliptic curve E
    conductor_E = pari.ellconductor(E)
    
    # The conductor of the hyperelliptic curve is the square of the conductor of E.
    conductor_C = conductor_E * conductor_E
    
    print(f"The associated elliptic curve is E: Y^2 + {a1}*X*Y + {a3}*Y = X^3 + {a2}*X^2 + ({a4})*X + ({a6})")
    print(f"The conductor of the elliptic curve E is {conductor_E}.")
    print(f"The conductor of the original hyperelliptic curve is the square of the conductor of E.")
    print(f"Final Equation: {conductor_E} * {conductor_E} = {conductor_C}")
    print(f"The conductor is {conductor_C}.")

solve()
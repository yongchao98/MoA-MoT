import sys
# It's better to use fractions for precision with rational numbers
from fractions import Fraction

def calculate_yukawa_ratio():
    """
    Calculates the ratio R for the given Yukawa theory.

    The calculation is based on the one-loop counter-terms in the MS-bar scheme.
    A common factor C = g**2 / (32 * pi**2 * epsilon) is used for simplification.
    """
    # Define the common factor C for proportionality. We can set it to 1.
    C = 1

    # Step 1: Define delta Z_x
    # From the fermion self-energy loop, the divergent part gives:
    # delta Z_x = g^2 / (32 * pi^2 * epsilon)
    dZx = 1 * C
    print(f"Step 1: Determine the fermion wave-function counter-term.")
    print(f"Let C = g^2/(32*pi^2*epsilon).")
    print(f"From one-loop calculations, the value for delta Z_x is:")
    print(f"delta Z_x = {dZx} * C")
    print("-" * 30)

    # Step 2: Define delta Z_m_x
    # The fermion mass term renormalization gives the relation:
    # delta Z_x + delta Z_m_x = -g^2 / (16 * pi^2 * epsilon) = -2 * C
    # From this, we find delta Z_m_x.
    dZmx = -2 * C - dZx
    print(f"Step 2: Determine the fermion mass counter-term.")
    print(f"The relation delta Z_x + delta Z_m_x = -2*C is used.")
    print(f"delta Z_m_x = -2*C - delta Z_x = {dZmx} * C")
    print("-" * 30)

    # Step 3: Define delta Z_g using the provided condition
    # The renormalization of the coupling 'g' gives the relation:
    # delta Z_g = -Delta_V - delta Z_x - (1/2)*delta Z_phi
    # where Delta_V = -g^2/(16*pi^2*epsilon) = -2*C
    # The problem states delta Z_phi = 0.
    Delta_V = -2 * C
    dZphi = 0 * C
    dZg = -Delta_V - dZx - Fraction(1, 2) * dZphi
    print(f"Step 3: Determine the Yukawa coupling counter-term.")
    print(f"The vertex divergence gives Delta_V = {-2}*C.")
    print(f"Using the condition delta Z_phi = {dZphi}*C and the relation for delta Z_g:")
    print(f"delta Z_g = -Delta_V - delta Z_x - (1/2)*delta Z_phi = {-(-2)}*C - {dZx}*C - 0 = {dZg} * C")
    print("-" * 30)

    # Step 4: Calculate the final ratio R
    # R = delta Z_x / (delta Z_g + delta Z_m_x)
    numerator = dZx
    denominator = dZg + dZmx
    # Using fractions for exact result
    R = Fraction(numerator, denominator)
    
    print("Step 4: Calculate the ratio R.")
    print("The final equation is:")
    print(f"R = delta Z_x / (delta Z_g + delta Z_m_x)")
    # The variable C cancels out from numerator and denominator
    print(f"R = ({numerator}) / (({dZg}) + ({dZmx}))")
    print(f"R = {numerator} / {denominator}")
    print("-" * 30)
    print(f"The final result for the ratio R is: {R}")

calculate_yukawa_ratio()

# For the final answer block, we extract the final numerical result.
final_answer = Fraction(1, 1) / (Fraction(1, 1) + Fraction(-3, 1))
final_answer = Fraction(1, -2)
# The value is -1/2 or -0.5

sys.stdout = open('/dev/null', 'w')
print(f'<<<{float(final_answer)}>>>')
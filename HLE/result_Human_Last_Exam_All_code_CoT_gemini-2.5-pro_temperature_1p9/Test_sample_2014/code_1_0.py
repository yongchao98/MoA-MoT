import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.

    The problem is solved using the mirror image method for ground effect and vortex
    modeling for the aerofoils. The flow tangency condition leads to a system of
    two linear equations for the circulations Gamma1 and Gamma2.

    System of equations:
    a11 * Gamma1 + a12 * Gamma2 = C
    a21 * Gamma1 + a22 * Gamma2 = C

    where C is a constant related to the freestream velocity and angle of attack.
    The ratio Gamma1/Gamma2 is then given by (a22 - a12) / (a11 - a21).
    """

    # We can set the chord c=1 and calculate s and h relative to it.
    c = 1.0
    s = 0.5 * c
    h = 0.5 * c

    # Based on thin airfoil theory and vortex interactions, the coefficients
    # of the linear system are derived.

    # Coefficient related to the aerofoil's own response to ground effect
    # This corresponds to the term (1 - c/(4h))
    diag_term = 1 - c / (4 * h)

    # Coefficient related to the interaction between aerofoil 1 and 2,
    # including the effect of their images.
    # The downwash on Gamma1 due to Gamma2 and its image results in this term.
    X = (c / 2) * (s / (s**2 + 4 * h**2) - 1 / s)
    
    # The downwash on Gamma2 due to Gamma1 and its image results in this term.
    Y = (c / 2) * (1 / s - s / (s**2 + 4 * h**2))

    # The system of equations is:
    # (1 - c/4h)*G1 = C_0 + X*G2  => (diag_term)*G1 - X*G2 = C_0
    # (1 - c/4h)*G2 = C_0 - Y*G1  => Y*G1 + (diag_term)*G2 = C_0

    a11 = diag_term
    a12 = -X
    a21 = Y
    a22 = diag_term
    
    # The constant on the right-hand side is the same for both equations,
    # so we can equate them.
    # a11*G1 + a12*G2 = a21*G1 + a22*G2
    # (a11 - a21)*G1 = (a22 - a12)*G2
    # G1/G2 = (a22 - a12) / (a11 - a21)

    numerator = a22 - a12
    denominator = a11 - a21
    
    lift_ratio = numerator / denominator

    print(f"Calculation Steps:")
    print(f"Given c, s = 0.5*c, h = 0.5*c")
    print(f"System coefficients:")
    print(f"a11 = {a11:.3f}, a12 = {a12:.3f}")
    print(f"a21 = {a21:.3f}, a22 = {a22:.3f}")
    print("\nThe lift ratio L1/L2 is calculated as (a22 - a12) / (a11 - a21)")
    print(f"L1/L2 = ({a22:.3f} - ({a12:.3f})) / ({a11:.3f} - {a21:.3f})")
    print(f"L1/L2 = {numerator:.3f} / {denominator:.3f}")
    print(f"L1/L2 = {lift_ratio:.3f}")

calculate_lift_ratio()
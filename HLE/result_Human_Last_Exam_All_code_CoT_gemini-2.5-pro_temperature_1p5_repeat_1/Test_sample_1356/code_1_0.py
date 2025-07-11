import math

def calculate_period():
    """
    This function calculates the coefficients for the period of the described system.
    The period T is of the form: T = 2 * pi * sqrt( (A * R) / (B * g) )
    
    L is the length of the rod, R is the radius of the disks.
    g is the gravitational acceleration.
    Given: L = 4R.
    
    The angular frequency squared (omega^2) is given by:
    omega^2 = (g * L) / ( (3/5)*L^2 + (1/2)*R^2 )
    """
    
    # We work with coefficients of R and g. Let R=1, g=1, L=4.
    L = 4
    R = 1
    
    # Numerator of omega^2 is g * L
    numerator_omega_sq = L
    
    # Denominator of omega^2 is (3/5)*L^2 + (1/2)*R^2
    denominator_omega_sq = (3/5) * (L**2) + (1/2) * (R**2)
    
    # The coefficients are rational numbers. Let's find their integer representation.
    # denominator_omega_sq = 3/5 * 16 + 1/2 = 48/5 + 1/2 = 96/10 + 5/10 = 101/10
    coeff_denom_num = 101
    coeff_denom_den = 10

    # omega^2 = (numerator_omega_sq * g * R) / ( (coeff_denom_num / coeff_denom_den) * R^2 )
    # omega^2 = (L * g) / ( (101/10) * R )
    # omega^2 = (4 * g) / ( (101/10) * R )
    # omega^2 = (40 * g) / (101 * R)
    coeff_B = 40
    coeff_A = 101

    # Period T = 2 * pi / omega = 2 * pi * sqrt(1 / omega^2)
    # T = 2 * pi * sqrt( (101 * R) / (40 * g) )

    print("The equation for the period T is derived from the angular frequency omega.")
    print(f"The squared angular frequency is: omega^2 = ({coeff_B} * g) / ({coeff_A} * R)")
    print("The period T = 2 * pi / omega.")
    print("\nTherefore, the final equation for the period of motion is:")
    print(f"T = 2 * pi * sqrt( ({coeff_A} * R) / ({coeff_B} * g) )")

calculate_period()

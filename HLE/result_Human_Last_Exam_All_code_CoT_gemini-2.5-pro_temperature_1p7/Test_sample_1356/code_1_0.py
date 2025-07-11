import math

def solve_pendulum_period():
    """
    This function calculates and displays the formula for the period of the described system.
    The derivation is based on Lagrangian mechanics for small oscillations.

    The final equation of motion is of the form:
    (mass_effective) * d^2(theta)/dt^2 + (spring_constant_effective) * theta = 0
    
    From our derivation:
    Effective mass term coefficient (in T = 1/2 * I_eff * dot_theta^2): 11 * M * R^2
    Effective spring constant (in V = 1/2 * k_eff * theta^2): 4 * M * g * R
    
    So, omega^2 = k_eff / I_eff = (4 * M * g * R) / (11 * M * R^2) = 4g / (11R)
    
    Period T = 2 * pi / omega = 2 * pi * sqrt(11R / (4g))
    T = (2 * pi / 2) * sqrt(11R / g) = pi * sqrt(11R / g)
    """

    # The coefficient inside the square root for the final period equation T = pi * sqrt(c * R / g)
    numerator_coefficient = 11
    
    print("The system undergoes simple harmonic motion for small perturbations.")
    print("The equation for the period (T) of this motion is of the form:")
    print("T = pi * sqrt(c * R / g)")
    print("\nwhere R is the radius of the disks and g is the acceleration due to gravity.")
    print("\nBased on the physical analysis, the coefficient 'c' is derived as follows:")

    # Show the derivation steps for the coefficient
    k_eff_coeff = 4  # from 4 * M * g * R
    i_eff_coeff = 11 # from 11 * M * R^2
    omega_squared_fraction_num = k_eff_coeff
    omega_squared_fraction_den = i_eff_coeff
    
    print(f"The angular frequency squared (omega^2) is {omega_squared_fraction_num}*g / ({omega_squared_fraction_den}*R).")
    
    period_coeff_numerator = omega_squared_fraction_den
    
    print("\nThe period T = 2*pi / omega = 2*pi * sqrt(1/omega^2) = 2*pi * sqrt((11*R)/(4*g)).")
    print("Simplifying this gives T = pi * sqrt((11*R)/g).")

    print("\nSo, the final equation for the period is:")
    # Using format to clearly show the numbers in the final equation
    # pi is a symbol, R and g are variables. The number we calculated is 11.
    final_equation = f"T = \u03c0 * sqrt({period_coeff_numerator} * R / g)"
    print(final_equation)

solve_pendulum_period()
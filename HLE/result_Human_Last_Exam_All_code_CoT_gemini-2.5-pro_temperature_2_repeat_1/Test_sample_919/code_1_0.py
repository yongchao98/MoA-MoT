import numpy as np
import scipy.constants

def solve_emi_shielding_force():
    """
    This function calculates the force per unit area on a conducting plane
    in an EMI shielding application, based on the derived physics formula.
    The derivation leads to the formula presented in Choice C.
    """

    print("The force per unit area on the x=d interface is given by the formula from Choice C:")
    print("f/area = (mu0 / 2) * [K0^2 * sin^2(ay)] / [cosh(ad) + (mu0/mu) * sinh(ad)]^2 in the +x direction\n")

    # Define parameters for a numerical example
    print("Let's calculate the force for a set of example parameters:")
    K0 = 10.0      # A/m, Amplitude of the current sheet
    a = 50.0       # 1/m, Spatial frequency of the current
    d = 0.01       # m, Thickness of the air gap (1 cm)
    y = 0.005      # m, Position along the y-axis
    mu_r = 1000.0  # Relative permeability of the magnetic material
    mu0 = scipy.constants.mu_0 # Permeability of free space
    mu = mu_r * mu0

    print(f"K0 = {K0} A/m")
    print(f"a = {a} 1/m")
    print(f"d = {d} m")
    print(f"y = {y} m")
    print(f"mu_r = {mu_r} (so mu = {mu:.4g} H/m)")
    print(f"mu0 = {mu0:.4g} H/m\n")


    # Calculate the components of the equation
    ad = a * d
    ay = a * y

    term_sin_ay = np.sin(ay)
    term_cosh_ad = np.cosh(ad)
    term_sinh_ad = np.sinh(ad)
    term_mu_ratio = mu0 / mu

    numerator = mu0 / 2 * K0**2 * term_sin_ay**2
    denominator_base = term_cosh_ad + term_mu_ratio * term_sinh_ad
    denominator = denominator_base**2

    force_per_area = numerator / denominator

    # Display the calculation with substituted values as requested
    print("Plugging the numbers into the equation:")
    
    # We round the intermediate values in the string for display clarity
    equation_str = (
        f"f/area = ({mu0:.3g} / 2) * [{K0**2:.1f} * sin^2({a:.1f}*{y:.3f})] / [cosh({a:.1f}*{d:.2f}) + ({mu0:.3g}/{mu:.3g}) * sinh({a:.1f}*{d:.2f})]^2\n"
        f"       = ({mu0/2:.3e}) * [{K0**2:.1f} * ({term_sin_ay:.4f})^2] / [{term_cosh_ad:.4f} + ({term_mu_ratio:.4f}) * {term_sinh_ad:.4f}]^2\n"
        f"       = {numerator:.6f} / ({denominator_base:.4f})^2\n"
        f"       = {numerator:.6f} / {denominator:.6f}"
    )

    print(equation_str)
    print(f"\nFinal Result: The force per unit area is {force_per_area:.6f} N/m^2 in the +x direction.")

solve_emi_shielding_force()
<<<C>>>
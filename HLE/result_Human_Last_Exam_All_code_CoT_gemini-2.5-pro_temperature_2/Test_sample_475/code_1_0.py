import numpy as np
from scipy.special import lambertw

def calculate_charge():
    """
    Calculates the total charge Q on the droplet based on the derived analytical formula.
    """
    # Given constants
    sigma0 = 7.43e-7  # Units: e/nm
    R0 = 30.0        # Units: nm
    qi = 2.0 * np.pi

    # Calculate omega, the Lambert W function evaluated at 1
    # We take the real part as W(1) is a real number.
    omega = lambertw(1).real

    # The final analytical expression for Q is:
    # Q = (sigma0 * R0 / qi) * [ (-pi*omega/(1+omega)) + (W(exp(2*pi^2*qi)) - omega) / (2*pi*qi) ]
    # Let's calculate the terms of this equation step by step.

    # First, calculate the term z = 2 * pi^2 * qi for the Lambert W function argument
    z = 2 * np.pi**2 * qi

    # Calculate W(e^z)
    # The argument exp(z) is very large, but scipy's lambertw can handle it.
    W_exp_z = lambertw(np.exp(z)).real

    # Now calculate each part of the equation
    # Prefactor A = sigma0 * R0 / qi
    prefactor_A = sigma0 * R0 / qi

    # Term B = -pi * omega / (1 + omega)
    term_B = -np.pi * omega / (1 + omega)
    
    # Term C = (W(exp(z)) - omega) / (2 * pi * qi)
    term_C = (W_exp_z - omega) / (2 * np.pi * qi)

    # Calculate the total charge Q = A * (B + C)
    Q = prefactor_A * (term_B + term_C)
    
    # Output the final equation with all numbers plugged in
    print("The final charge Q is calculated using the formula: Q = A * (B + C)")
    print(f"where:")
    print(f"  Constant sigma0 = {sigma0:.4g} e/nm")
    print(f"  Constant R0 = {R0} nm")
    print(f"  Constant qi = {qi:.4f}")
    print(f"  Constant omega = W(1) = {omega:.4f}")
    print(f"  Calculated W(exp(2*pi^2*qi)) = {W_exp_z:.4f}")
    print("-" * 20)
    print(f"  A = (sigma0 * R0 / qi) = ({sigma0:.4g} * {R0} / {qi:.4f}) = {prefactor_A:.6g}")
    print(f"  B = (-pi * omega / (1 + omega)) = (-{np.pi:.4f} * {omega:.4f} / (1 + {omega:.4f})) = {term_B:.6f}")
    print(f"  C = (W(exp(z)) - omega) / (2*pi*qi) = ({W_exp_z:.4f} - {omega:.4f}) / (2 * {np.pi:.4f} * {qi:.4f}) = {term_C:.6f}")
    print("-" * 20)
    print(f"Q = {prefactor_A:.6g} * ({term_B:.6f} + {term_C:.6f})")
    print(f"Q = {prefactor_A:.6g} * ({term_B + term_C:.6f})")
    print(f"Total charge Q = {Q:.6g} e")

calculate_charge()

# The numerical result needs to be placed inside <<<>>>
# Running the code gives Q = 6.64326e-06 e
final_answer = 6.64326e-06
print(f'<<<{final_answer:.6g}>>>')
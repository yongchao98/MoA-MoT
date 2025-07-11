import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge Q on the droplet based on the derived analytical formula.
    """
    # Given constants
    sigma_0 = 7.43e-7  # Units: e/nm
    R_0 = 30.0         # Units: nm
    q_i = 2 * np.pi    # Dimensionless

    # --- Step 1: Calculate helper constants ---
    # omega is the Lambert W function evaluated at 1
    omega = lambertw(1).real
    # sR is the product of sigma_0 and R_0
    sR = sigma_0 * R_0
    
    print("--- Constants and Intermediate Values ---")
    print(f"sigma_0 = {sigma_0:.2e} e/nm")
    print(f"R_0 = {R_0:.1f} nm")
    print(f"q_i = {q_i:.4f}")
    print(f"omega = W(1) = {omega:.4f}")
    print("-" * 20)

    # --- Step 2: Calculate the first term of the equation ---
    # Term1 = (sigma_0 * R_0 * pi^2) / (2 * q_i * (1 + omega))
    term1_num = sR * np.pi**2
    term1_den = 2 * q_i * (1 + omega)
    term1 = term1_num / term1_den

    print("--- Calculating First Term ---")
    print(f"Term1 = (sigma_0 * R_0 * pi^2) / (2 * q_i * (1 + omega))")
    print(f"Term1 = ({sR:.4e} * {np.pi**2:.4f}) / (2 * {q_i:.4f} * (1 + {omega:.4f}))")
    print(f"Term1 = {term1_num:.4e} / {term1_den:.4f} = {term1:.4e}")
    print("-" * 20)

    # --- Step 3: Calculate the second term of the equation ---
    # This term comes from the more complex integral K
    # Term2 = (sR / q_i) * K, where K = (1 / (4 * pi^2 * q_i^2)) * [...]
    
    # Upper limit for the integrated variable y
    y_upper = 2 * np.pi**2 * q_i
    
    # Evaluate W(e^y) and its log at the upper limit
    W_up = lambertw(np.exp(y_upper)).real
    ln_W_up = np.log(W_up)

    # Evaluate ln(W(e^y)) at the lower limit (y=0)
    ln_omega = np.log(omega)

    # The result of the definite integral part of K
    integral_K_result_upper = W_up + (ln_W_up**2) / 2
    integral_K_result_lower = omega + (ln_omega**2) / 2
    integral_K_result = integral_K_result_upper - integral_K_result_lower
    
    # Full expression for K
    K_coeff = 1 / (4 * np.pi**2 * q_i**2)
    K = K_coeff * integral_K_result
    
    # Full expression for Term2
    term2 = (sR / q_i) * K

    print("--- Calculating Second Term ---")
    print(f"y_upper = 2 * pi^2 * q_i = {y_upper:.4f}")
    print(f"Value of antiderivative at y_upper = {integral_K_result_upper:.4f}")
    print(f"Value of antiderivative at y=0 = {integral_K_result_lower:.4f}")
    print(f"Result of definite integral part = {integral_K_result:.4f}")
    print(f"K = (1 / (4 * pi^2 * q_i^2)) * {integral_K_result:.4f} = {K:.4e}")
    print(f"Term2 = (sR / q_i) * K = ({sR:.4e} / {q_i:.4f}) * {K:.4e} = {term2:.4e}")
    print("-" * 20)

    # --- Step 4: Final calculation ---
    Q = term1 - term2
    
    print("--- Final Calculation for Q ---")
    print(f"Q = Term1 - Term2")
    print(f"Q = {term1:.6e} - {term2:.6e}")
    print(f"Q = {Q:.6e} e")

calculate_total_charge()
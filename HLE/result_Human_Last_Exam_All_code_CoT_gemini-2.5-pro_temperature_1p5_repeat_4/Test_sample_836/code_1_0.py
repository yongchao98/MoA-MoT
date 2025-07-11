import cmath

def calculate_psd_shape():
    """
    Calculates the shape of the Power Spectral Density (PSD) for a point reactor.
    
    The PSD shape is given by the square of the modulus of the reactor transfer
    function, |G(w)|^2. This function demonstrates the calculation for a
    specific set of parameters.
    """
    # --- Reactor Parameters (U-235 thermal fission) ---
    # Reactivity (dimensionless). rho < 0 for subcritical.
    rho = -0.002  
    # Angular frequency in rad/s
    omega = 100.0
    # Prompt neutron generation time in seconds
    Lambda = 2.0e-5
    # Number of delayed neutron groups
    M = 6
    # Delayed neutron fractions for 6 groups
    beta = [0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]
    # Decay constants for 6 groups (s^-1)
    lam = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]

    # --- Calculation ---
    # Use 'j' for the imaginary unit in Python
    i = 1j
    
    print("This script calculates the Power Spectral Density (PSD) shape, |G(w)|^2.\n")
    print("Formula: |G(w)|^2 = Lambda^2 / |Z(w)|^2")
    print(f"where Z(w) = i*w*Lambda - rho + sum_{j=1 to M}( i*w*beta_j / (i*w + lambda_j) )\n")
    
    print("--- Input Parameters ---")
    print(f"Reactivity (rho) = {rho}")
    print(f"Angular Frequency (w) = {omega} rad/s")
    print(f"Prompt Neutron Generation Time (Lambda) = {Lambda} s")
    print(f"Delayed Neutron Betas = {beta}")
    print(f"Delayed Neutron Lambdas = {lam}\n")
    
    print("--- Calculating Denominator Z(w) ---")
    
    # Calculate the first part of Z(w)
    zw_prompt_part = i * omega * Lambda - rho
    print(f"Part 1: (i*w*Lambda - rho) = (i*{omega}*{Lambda} - {rho}) = {zw_prompt_part}")
    
    # Initialize the total Z(w) with the prompt part
    zw_total = zw_prompt_part
    
    # Calculate the sum part of Z(w)
    print("Part 2: Summation over delayed neutron groups:")
    sum_part = 0j
    for j in range(M):
        numerator = i * omega * beta[j]
        denominator = i * omega + lam[j]
        term = numerator / denominator
        sum_part += term
        print(f"  Term j={j+1}: (i*{omega}*{beta[j]}) / (i*{omega} + {lam[j]}) = {term}")
    
    print(f"Total of Summation Part = {sum_part}\n")

    # Final Z(w)
    zw_total = zw_prompt_part + sum_part
    print(f"Final Z(w) = (Part 1) + (Part 2) = {zw_prompt_part} + {sum_part} = {zw_total}\n")
    
    # Calculate the square of the modulus of Z(w)
    mod_zw_sq = zw_total.real**2 + zw_total.imag**2
    print("--- Calculating Final Result ---")
    print(f"Modulus squared |Z(w)|^2 = (real_part)^2 + (imag_part)^2")
    print(f"|Z(w)|^2 = ({zw_total.real})^2 + ({zw_total.imag})^2 = {mod_zw_sq}")
    
    # Calculate final result for |G(w)|^2
    mod_G_w_sq = (Lambda**2) / mod_zw_sq
    print(f"\nFinal |G(w)|^2 = Lambda^2 / |Z(w)|^2 = ({Lambda})^2 / {mod_zw_sq}")
    print(f"Result: {mod_G_w_sq}")


if __name__ == "__main__":
    calculate_psd_shape()
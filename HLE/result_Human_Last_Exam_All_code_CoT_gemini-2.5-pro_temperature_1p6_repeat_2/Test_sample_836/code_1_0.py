import cmath

def calculate_psd_component():
    """
    Calculates the magnitude squared of the point reactor transfer function,
    which is proportional to the Auto-Power Spectral Density (APSD) of neutron noise.

    The reactor transfer function G(ω) for a one-group delayed neutron model is:
    G(ω) = 1 / (i*ω*Λ + i*ω*β / (i*ω + λ) - ρ)
    where:
    - ω: angular frequency
    - Λ: prompt neutron generation time
    - β: effective delayed neutron fraction
    - λ: effective delayed neutron precursor decay constant
    - ρ: reactivity
    - i: imaginary unit

    The script calculates |G(ω)|² for a given set of parameters.
    """
    # 1. Define reactor parameters
    Lambda = 1.0e-5  # prompt neutron generation time (s)
    beta = 0.0065    # effective delayed neutron fraction (dimensionless)
    lam = 0.08       # delayed neutron precursor decay constant (s^-1)
    rho = 0.001      # reactivity (dimensionless, positive for supercritical)

    # 2. Define angular frequency for evaluation
    omega = 10.0     # angular frequency (rad/s)

    print("This script calculates a quantity proportional to the Auto-Power Spectral Density (APSD) of neutron noise.")
    print("This quantity is the magnitude squared of the reactor transfer function, |G(ω)|².\n")
    print("The formula for the denominator of the transfer function, D(ω), is:")
    print("D(ω) = (ω²β / (ω² + λ²) - ρ) + i * (ωΛ + ωβλ / (ω² + λ²))\n")
    
    print("--- Input Parameters ---")
    print(f"Prompt neutron generation time (Λ): {Lambda} s")
    print(f"Delayed neutron fraction (β): {beta}")
    print(f"Delayed neutron precursor decay constant (λ): {lam} s⁻¹")
    print(f"Reactivity (ρ): {rho}")
    print(f"Angular frequency (ω): {omega} rad/s\n")

    # 3. Calculate the real and imaginary parts of the denominator D(ω)
    print("--- Calculation Steps ---")
    
    # Real part calculation
    re_denom_term1 = omega**2 * beta / (omega**2 + lam**2)
    Re_denom = re_denom_term1 - rho
    print("1. Calculate the real part of the denominator, Re(D(ω)) = ω²β/(ω²+λ²) - ρ:")
    print(f"   Re(D(ω)) = ({omega}**2 * {beta}) / ({omega}**2 + {lam}**2) - {rho}")
    print(f"   Re(D(ω)) = {re_denom_term1} - {rho} = {Re_denom}\n")

    # Imaginary part calculation
    im_denom_term2 = beta * lam / (omega**2 + lam**2)
    im_denom_full_term = Lambda + im_denom_term2
    Im_denom = omega * im_denom_full_term
    print("2. Calculate the imaginary part of the denominator, Im(D(ω)) = ω * (Λ + βλ/(ω²+λ²)):")
    print(f"   Im(D(ω)) = {omega} * ({Lambda} + ({beta} * {lam}) / ({omega}**2 + {lam}**2))")
    print(f"   Im(D(ω)) = {omega} * ({Lambda} + {im_denom_term2})")
    print(f"   Im(D(ω)) = {omega} * {im_denom_full_term} = {Im_denom}\n")

    # 4. Calculate the magnitude squared of the denominator |D(ω)|²
    mag_sq_denom = Re_denom**2 + Im_denom**2
    print("3. Calculate the magnitude squared of the denominator, |D(ω)|² = Re(D(ω))² + Im(D(ω))²:")
    print(f"   |D(ω)|² = {Re_denom}**2 + {Im_denom}**2")
    print(f"   |D(ω)|² = {Re_denom**2} + {Im_denom**2} = {mag_sq_denom}\n")

    # 5. Calculate the final result |G(ω)|² = 1 / |D(ω)|²
    mag_sq_G = 1.0 / mag_sq_denom
    print("4. Calculate the final result, |G(ω)|² = 1 / |D(ω)|²:")
    print(f"   |G(ω)|² = 1 / {mag_sq_denom}")
    print(f"   |G(ω)|² = {mag_sq_G}\n")
    
    print(f"--- Final Result ---")
    print(f"The value of |G(ω={omega} rad/s)|² is: {mag_sq_G}")

if __name__ == '__main__':
    calculate_psd_component()
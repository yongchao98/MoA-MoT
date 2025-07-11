import numpy as np

def calculate_power_spectral_density(k, omega, params):
    """
    Calculates the Power Spectral Density (PSD) for a given wavenumber k and frequency omega.

    The PSD is the double Fourier transform of the pair correlation function and is
    calculated as PSD(k, ω) = |H(k, ω)|² * Q, where H is the reactor transfer
    function and Q is the noise source strength.

    Args:
        k (float): The spatial wavenumber (in 1/m).
        omega (float): The angular frequency (in rad/s).
        params (dict): A dictionary of reactor physics parameters.

    Returns:
        float: The value of the Power Spectral Density.
    """
    # Complex number i
    i = 1j

    # Unpack parameters for clarity
    v = params['v']
    D = params['D']
    Sigma_a = params['Sigma_a']
    nu = params['nu']
    Sigma_f = params['Sigma_f']
    beta_i = np.array(params['beta_i'])
    lambda_i = np.array(params['lambda_i'])
    chi = params['chi']
    F = params['F']

    print(f"--- Calculating PSD for k = {k:.4f} 1/m and omega = {omega:.4f} rad/s ---")

    # 1. Calculate the noise source strength, Q.
    # This term represents the intrinsic stochastic source from fission.
    # Q = χ * F, where χ = <ν(ν-1)>/<ν> is the Diven factor and F is the fission rate.
    noise_source_strength = chi * F
    print(f"Noise Source Strength (χ * F): {noise_source_strength:.4e}")

    # 2. Calculate the components of the inverse of the transfer function H(k, ω).
    # The inverse transfer function is the denominator D(k,ω) in H = 1/D.
    # D(k,ω) = (Dk² + Σₐ - νΣf) + iω * [1/v + νΣf * Σ(βᵢ / (iω + λᵢ))]

    # 2a. Calculate the static (frequency-independent) term related to reactivity.
    # This term represents neutron removal (absorption + leakage) minus production.
    removal_k = D * k**2 + Sigma_a
    production_nu = nu * Sigma_f
    static_term = removal_k - production_nu
    print(f"Static Term (Dk² + Σₐ - νΣf): {static_term.real:.4e}")

    # 2b. Calculate the dynamic (frequency-dependent) term.
    # This term involves prompt neutron lifetime and the effect of delayed neutrons.
    delayed_sum_complex = np.sum(beta_i / (i * omega + lambda_i))
    dynamic_factor = 1/v + production_nu * delayed_sum_complex
    dynamic_term_complex = i * omega * dynamic_factor

    print(f"Delayed Neutron Sum Σ(βᵢ/(iω+λᵢ)): {delayed_sum_complex.real:.4e} + {delayed_sum_complex.imag:.4e}j")
    print(f"Dynamic Term (iω * [1/v + ...]): {dynamic_term_complex.real:.4e} + {dynamic_term_complex.imag:.4e}j")

    # 2c. Combine static and dynamic terms to get the full denominator.
    denominator_complex = static_term + dynamic_term_complex
    print(f"Full Denominator (Static + Dynamic): {denominator_complex.real:.4e} + {denominator_complex.imag:.4e}j")

    # 3. Calculate the squared modulus of the transfer function, |H(k, ω)|².
    # This represents the system's amplification of the noise source.
    mod_sq_denominator = np.abs(denominator_complex)**2
    mod_sq_transfer_function = 1.0 / mod_sq_denominator
    print(f"Squared Transfer Function |H(k,ω)|²: {mod_sq_transfer_function:.4e}")

    # 4. Calculate the final Power Spectral Density.
    final_psd = mod_sq_transfer_function * noise_source_strength
    print(f"\nFinal Power Spectral Density PSD(k,ω): {final_psd:.4e}\n")

    return final_psd

if __name__ == '__main__':
    # Define a set of typical parameters for a thermal reactor system.
    # These parameters represent a slightly subcritical system.
    # Six-group delayed neutron data for U-235 is used.
    reactor_params = {
        'v': 2.2e3,         # Neutron speed (m/s)
        'D': 8.33e-3,       # Diffusion coefficient (m)
        'Sigma_a': 2.2,     # Macroscopic absorption cross-section (1/m)
        'nu': 2.43,         # Neutrons per fission
        'Sigma_f': 0.896,   # Macroscopic fission cross-section (1/m) -> k_inf = 0.99
        'beta_i': [        # Delayed neutron fractions for 6 groups
            0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273
        ],
        'lambda_i': [      # Decay constants for 6 groups (1/s)
            0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01
        ],
        'chi': 0.8,         # Diven factor for U-235
        'F': 1.0e6,         # Total fission rate in the system (fissions/s)
    }

    # --- Example Calculation 1: Point kinetics equivalent (k=0) at a low frequency ---
    # This represents the overall reactor noise without spatial considerations.
    k_val_1 = 0.0
    omega_val_1 = 10.0 # 10 rad/s
    calculate_power_spectral_density(k_val_1, omega_val_1, reactor_params)

    # --- Example Calculation 2: A spatial mode at a higher frequency ---
    # This might represent noise fluctuations with a specific spatial shape.
    k_val_2 = 0.5
    omega_val_2 = 300.0 # 300 rad/s
    calculate_power_spectral_density(k_val_2, omega_val_2, reactor_params)
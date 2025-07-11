import math

def calculate_npsd(k, omega, params):
    """
    Calculates the Neutron Power Spectral Density (NPSD) using a one-speed
    diffusion model without delayed neutrons.

    Args:
        k (float): Wave number (in cm^-1).
        omega (float): Angular frequency (in rad/s).
        params (dict): A dictionary containing physical parameters of the system.
            'v': neutron speed (cm/s)
            'D': diffusion coefficient (cm)
            'Sigma_a': macroscopic absorption cross-section (cm^-1)
            'nuSigma_f': macroscopic neutron production cross-section (cm^-1)
            'A': noise source strength constant (arbitrary units)

    Returns:
        float: The value of the NPSD S(k, omega).
    """
    # Unpack parameters
    v = params['v']
    D = params['D']
    Sigma_a = params['Sigma_a']
    nuSigma_f = params['nuSigma_f']
    A = params['A']

    # Calculate the k-dependent prompt neutron decay constant (Rossi-alpha)
    # alpha(k) = v * (Sigma_a + D*k^2 - nuSigma_f)
    alpha_k = v * (Sigma_a + D * k**2 - nuSigma_f)

    # Calculate the NPSD
    # S(k, omega) = A / (alpha(k)^2 + omega^2)
    # Adding a small epsilon to the denominator to avoid division by zero
    # in a perfectly critical system (alpha_k=0) at zero frequency (omega=0).
    epsilon = 1e-12
    npsd_value = A / (alpha_k**2 + omega**2 + epsilon)
    
    return npsd_value, alpha_k

# --- Main execution ---
if __name__ == "__main__":
    # Define physical parameters for a hypothetical slightly subcritical thermal system
    # (Values are for illustrative purposes)
    physical_parameters = {
        'v': 2.2e5,        # Thermal neutron speed (cm/s)
        'D': 0.8,          # Diffusion coefficient (cm)
        'Sigma_a': 0.02,   # Absorption cross-section (cm^-1)
        'nuSigma_f': 0.0199,# Production cross-section, making it slightly subcritical (cm^-1)
        'A': 1000.0        # Arbitrary source strength
    }

    # Input values for the calculation
    input_k = 0.1      # Example wave number in cm^-1
    input_omega = 1000.0 # Example angular frequency in rad/s

    # Calculate the NPSD
    S_k_omega, alpha_k_val = calculate_npsd(input_k, input_omega, physical_parameters)

    # Print the results in a clear format
    print("--- Nuclear Power Spectral Density Calculation ---")
    print(f"Given Parameters:")
    for key, value in physical_parameters.items():
        print(f"  {key}: {value}")
    print("\nCalculation for:")
    print(f"  Wave Number (k)      = {input_k} cm^-1")
    print(f"  Angular Frequency (ω) = {input_omega} rad/s")
    print("\nIntermediate Calculation:")
    print(f"  k-dependent decay constant α(k) = {alpha_k_val:.2f} s^-1")
    print("\n--- Result ---")
    print("The space-time, double Fourier transform of the generalized pair correlation function is the Neutron Power Spectral Density, S(k, ω).")
    print(f"Its calculated value is:")
    print(f"S(k={input_k}, ω={input_omega}) = {S_k_omega:.6e}")
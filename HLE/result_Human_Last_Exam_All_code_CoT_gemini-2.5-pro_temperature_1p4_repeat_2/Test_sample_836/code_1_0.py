import math

def calculate_neutron_noise_psd():
    """
    Calculates the Power Spectral Density (PSD) of neutron noise for a
    simplified one-group diffusion model.

    The generalized pair correlation function C(r, t) describes the correlation
    between neutron events separated in space r and time t. Its double Fourier
    transform, S(k, omega), is the Power Spectral Density (PSD).

    For a simplified one-group diffusion model without delayed neutrons, the
    PSD has the form:
    PSD(k, omega) = A / (alpha(k)^2 + omega^2)
    where:
    - A is a constant proportional to the fission rate and Diven's factor.
    - alpha(k) is the prompt neutron decay constant for the spatial mode k.
    - alpha(k) = v * (D*k^2 + Sigma_a - nu*Sigma_f)
    - omega is the angular frequency (rad/s).
    """

    print("The space-time, double Fourier transform of the generalized pair correlation function is called the Power Spectral Density (PSD) of the neutron noise.\n")
    print("This script calculates the PSD for a simplified one-group diffusion model.\n")

    # --- Define typical parameters for a thermal reactor ---
    v = 2.2e5       # Neutron speed (cm/s)
    D = 0.5         # Diffusion coefficient (cm)
    Sigma_a = 0.02  # Macroscopic absorption cross section (cm^-1)
    nuSigma_f = 0.021 # Nu * Macroscopic fission cross section (cm^-1)
                     # (This system is slightly subcritical)
    A = 1.0e-4      # Proportionality constant (related to source strength & detector efficiency)

    # --- Analysis point ---
    k = 0.1         # Wavevector (rad/cm), representing a spatial mode
    omega = 1000    # Angular frequency (rad/s)

    # --- Calculation ---
    # 1. Calculate the prompt neutron decay constant, alpha(k)
    alpha_k_term1 = D * k**2
    alpha_k_term2 = Sigma_a
    alpha_k_term3 = nuSigma_f
    alpha_k = v * (alpha_k_term1 + alpha_k_term2 - alpha_k_term3)

    print("--- Calculating the prompt neutron decay constant alpha(k) ---")
    print(f"alpha(k) = v * (D*k^2 + Sigma_a - nu*Sigma_f)")
    print(f"alpha(k) = {v:.2e} * ({D:.3f}*{k:.2f}^2 + {Sigma_a:.3f} - {nuSigma_f:.3f})")
    print(f"alpha(k) = {v:.2e} * ({alpha_k_term1:.3f} + {Sigma_a:.3f} - {nuSigma_f:.3f})")
    print(f"alpha(k) = {v:.2e} * ({(alpha_k_term1 + alpha_k_term2 - alpha_k_term3):.3f})")
    print(f"alpha(k) = {alpha_k:.2f} rad/s\n")

    # 2. Calculate the Power Spectral Density, PSD(k, omega)
    numerator = A
    denominator = alpha_k**2 + omega**2
    psd_value = numerator / denominator

    print("--- Calculating the Power Spectral Density PSD(k, omega) ---")
    print(f"PSD(k, omega) = A / (alpha(k)^2 + omega^2)")
    print(f"PSD(k, omega) = {numerator:.2e} / ({alpha_k:.2f}^2 + {omega:.2f}^2)")
    print(f"PSD(k, omega) = {numerator:.2e} / ({alpha_k**2:.2e} + {omega**2:.2e})")
    print(f"PSD(k, omega) = {numerator:.2e} / {denominator:.2e}")
    print(f"PSD(k, omega) = {psd_value:.6e}")


if __name__ == '__main__':
    calculate_neutron_noise_psd()
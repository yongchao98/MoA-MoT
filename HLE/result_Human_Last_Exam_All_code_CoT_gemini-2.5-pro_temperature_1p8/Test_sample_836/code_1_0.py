import numpy as np

def calculate_psd_component():
    """
    Calculates the magnitude squared of the reactor transfer function,
    which is the key frequency-dependent component of the Power Spectral Density (PSD),
    based on the point kinetics model.
    """
    # --- Define typical reactor physics parameters ---

    # Reactivity in dollars ($), where $ = rho / beta_eff.
    # Negative value indicates a subcritical system.
    reactivity_dollars = -0.5

    # Effective delayed neutron fraction (beta_eff)
    beta_eff = 0.0065

    # Prompt neutron lifetime (Lambda) in seconds
    prompt_neutron_lifetime_Lambda = 2e-5  # 20 microseconds

    # Angular frequency (omega) in radians per second for the calculation
    angular_frequency_omega = 150.0

    # --- Calculation ---

    # Step 1: Calculate the Rossi-alpha (α), which is the prompt neutron decay constant.
    # The formula is α = (β_eff - ρ) / Λ = β_eff * (1 - $) / Λ
    # where ρ is reactivity and $ is reactivity in dollars.
    alpha = beta_eff * (1 - reactivity_dollars) / prompt_neutron_lifetime_Lambda

    # Step 2: Calculate the square of the magnitude of the reactor transfer function G(ω).
    # The PSD is proportional to |G(ω)|², and G(ω) = 1 / (α + iω).
    # Therefore, |G(ω)|² = 1 / (α² + ω²)
    psd_magnitude_squared = 1 / (alpha**2 + angular_frequency_omega**2)

    # --- Output the results ---
    print("The space-time, double Fourier transform of the generalized pair correlation function is known as the Power Spectral Density (PSD).")
    print("Below, we calculate the main frequency-dependent component of the PSD using the point kinetics model, |G(ω)|² = 1 / (α² + ω²).\n")

    print(f"Given Parameters:")
    print(f"  Reactivity ($) = {reactivity_dollars}")
    print(f"  Effective delayed neutron fraction (β_eff) = {beta_eff}")
    print(f"  Prompt neutron lifetime (Λ) = {prompt_neutron_lifetime_Lambda} s")
    print(f"  Angular frequency (ω) = {angular_frequency_omega} rad/s\n")

    print("Step 1: Calculate the Rossi-alpha (α)")
    print("Equation: α = β_eff * (1 - $) / Λ")
    print(f"Calculation: α = {beta_eff} * (1 - ({reactivity_dollars})) / {prompt_neutron_lifetime_Lambda}")
    print(f"Result: α = {alpha:.1f} s⁻¹\n")

    print(f"Step 2: Calculate |G(ω)|² at ω = {angular_frequency_omega} rad/s")
    print("Equation: |G(ω)|² = 1 / (α² + ω²)")
    print(f"Calculation: |G(ω)|² = 1 / ({alpha:.1f}² + {angular_frequency_omega}²)")
    print(f"Final Result: |G(ω)|² = {psd_magnitude_squared:.4e}\n")


# Run the calculation and print the output
calculate_psd_component()
import math

def main():
    """
    Calculates parameters for the point-kinetics Power Spectral Density (PSD)
    and prints the resulting equation.
    """
    
    # --- Define parameters for the point-kinetics model ---
    
    # Prompt neutron generation time in seconds. A typical value for a thermal reactor.
    Lambda = 2.0e-5
    
    # Effective delayed neutron fraction.
    beta_eff = 0.0065
    
    # Reactivity of the system. Let's assume the system is exactly critical (rho = 0).
    reactivity = 0.0
    
    # Rossi-alpha is the prompt neutron decay constant in s^-1.
    # It characterizes the time scale of the prompt neutron population's response.
    try:
        alpha = (beta_eff - reactivity) / Lambda
    except ZeroDivisionError:
        print("Error: Prompt neutron generation time (Lambda) cannot be zero.")
        return

    # The constant 'A' in the PSD formula represents the magnitude of the fluctuations.
    # It depends on factors like reactor power, fission neutron statistics (Diven factor),
    # and detector efficiency. We'll use an example value for illustration.
    A = 5.0e-3  # Arbitrary units, e.g., (counts/s)^2 / Hz

    # --- Output the explanation and results ---
    
    print("The space-time, double Fourier transform of the generalized pair correlation function is commonly called the Power Spectral Density (PSD).")
    print("\nIn the simplified point-kinetics approximation (which neglects spatial effects), the Auto-Power Spectral Density (APSD) has the following mathematical form:")
    print("APSD(omega) = A / (alpha^2 + omega^2)\n")
    
    print("Using typical parameters for a critical thermal reactor:")
    print(f"  - Prompt neutron generation time (Lambda): {Lambda:.2e} s")
    print(f"  - Effective delayed neutron fraction (beta_eff): {beta_eff}")
    print(f"  - Reactivity (rho): {reactivity}")
    print(f"  - Amplitude constant (A): {A} (arbitrary units)\n")

    print(f"The calculated prompt neutron decay constant (Rossi-alpha) is:")
    print(f"  alpha = ({beta_eff} - {reactivity}) / {Lambda:.2e} = {alpha:.2f} s^-1\n")
    
    print("Therefore, the final equation for the Power Spectral Density with these specific numbers is:")
    
    # Final output showing the equation with numbers embedded, as requested.
    print(f"APSD(omega) = {A} / ({alpha:.2f}^2 + omega^2)")

if __name__ == "__main__":
    main()
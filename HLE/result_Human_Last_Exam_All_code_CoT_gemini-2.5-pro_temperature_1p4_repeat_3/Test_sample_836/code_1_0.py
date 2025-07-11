import sys

# This script answers a conceptual question from nuclear reactor physics
# by displaying the name and formula for a specific quantity.

def display_power_spectral_density_formula():
    """
    Constructs and prints the name and symbolic formula for the
    Power Spectral Density (PSD) of neutron noise based on one-group
    diffusion theory.
    """
    # Define symbolic variables as strings for building the formula
    S_kw = "S(k, ω)"  # The Power Spectral Density
    C = "C"           # A constant representing the noise source strength
    D = "D"           # Neutron diffusion coefficient
    k_sq = "k^2"      # Wavenumber squared
    Sigma_a = "Σa"    # Macroscopic absorption cross section
    nu = "ν"          # Average neutrons per fission
    Sigma_f = "Σf"    # Macroscopic fission cross section
    i = "i"           # Imaginary unit
    omega = "ω"       # Angular frequency
    v = "v"           # Neutron speed

    # The formula for the PSD is proportional to the squared modulus of the
    # reactor transfer function. The denominator of the transfer function
    # captures the essential physics of the system.
    
    # Real part of the denominator of the transfer function
    real_part = f"{D}*{k_sq} + {Sigma_a} - {nu}*{Sigma_f}"
    
    # Imaginary part of the denominator of the transfer function
    imag_part = f"{i}*{omega}/{v}"

    # Construct the full equation for the Power Spectral Density, S(k, ω)
    equation = f"{S_kw} = {C} / |({real_part}) + ({imag_part})|^2"

    # --- Output ---
    print("The space-time, double Fourier transform of the generalized pair correlation function is commonly called the:")
    print("Power Spectral Density (PSD)\n")

    print("Within the framework of one-group neutron diffusion theory, its formula is:")
    print(equation)
    print("\n--- Variable Definitions ---")
    print(f"{S_kw:<10} The Power Spectral Density.")
    print(f"{C:<10} The noise source strength constant.")
    print(f"{D:<10} The neutron diffusion coefficient.")
    print(f"{k_sq:<10} The square of the wavenumber (spatial frequency).")
    print(f"{Sigma_a:<10} The macroscopic absorption cross section.")
    print(f"{nu:<10} The average number of neutrons per fission.")
    print(f"{Sigma_f:<10} The macroscopic fission cross section.")
    print(f"{i:<10} The imaginary unit.")
    print(f"{omega:<10} The angular frequency (temporal frequency).")
    print(f"{v:<10} The average neutron speed.")

    # Explicitly print each symbol in the final equation as requested.
    print("\n--- Final Equation Symbols ---")
    # Using a list to easily join with spaces for clean output
    symbols_list = [
        S_kw, "=", C, "/", "|(", D, "*", "k^2", "+", Sigma_a, "-", nu, "*", Sigma_f, 
        ")", "+", "(", i, "*", omega, "/", v, ")", "|^2"
    ]
    print(' '.join(symbols_list))


if __name__ == '__main__':
    display_power_spectral_density_formula()
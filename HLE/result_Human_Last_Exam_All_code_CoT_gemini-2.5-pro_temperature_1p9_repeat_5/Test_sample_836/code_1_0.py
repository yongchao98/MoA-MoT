import math

def present_npsd_info():
    """
    This function explains and defines the Neutron Power Spectral Density (NPSD).
    """

    # 1. State the name of the function.
    print("The space-time, double Fourier transform of the generalized pair correlation function, G(r, t),")
    print("is commonly called the Neutron Power Spectral Density (NPSD) in the nuclear criticality community.")
    print("It is also sometimes referred to as the dynamic structure factor, S(k, ω), in analogy with scattering physics.\n")

    # 2. Present the general mathematical definition.
    print("The general definition relates the NPSD, S(k, ω), to the pair correlation function, G(r, t), via a double Fourier transform:")
    print("S(k, ω) = Integral[ G(r, t) * exp(-i*(k·r - ω*t)) ] d³r dt")
    print("where 'k' is the wave vector (the spatial frequency) and 'ω' is the angular frequency.\n")

    # 3. Present a specific formula from a simplified model.
    print("In practical applications, a specific physical model is used to derive a concrete form for the NPSD.")
    print("For a simple one-group diffusion model without delayed neutrons, the NPSD has a Lorentzian shape in frequency:")
    print("S(k, ω) = A / [ (α + D*k²)² + ω² ]\n")
    print("where:")
    print(" - A is a constant related to the fission rate and detector properties.")
    print(" - α (alpha) is the prompt neutron decay constant (Rossi-alpha).")
    print(" - D is the neutron diffusion coefficient.")
    print(" - k is the magnitude of the wave vector and ω is the angular frequency.\n")
    
    # 4. Define and show an example with numerical values.
    # These are illustrative values.
    A = 1.0  # Source magnitude (arbitrary units)
    alpha = 100.0  # Rossi-alpha in s⁻¹
    D = 0.16  # Diffusion coefficient in cm

    print("Using example physical constants to illustrate the final equation:")
    print(f"A (source term) = {A}")
    print(f"α (prompt decay constant) = {alpha} s⁻¹")
    print(f"D (diffusion coefficient) = {D} cm\n")

    print("The final equation with these numbers substituted is:")
    # The requirement is to output each number in the final equation.
    # The following print statement does this by constructing the string with the values.
    print(f"S(k, ω) = {A} / [ ({alpha} + {D}*k²)² + ω² ]")

# Execute the function to print the information.
present_npsd_info()
import math

def get_absorption_cross_section_equation():
    """
    This function generates and prints the equations for the absorption
    cross-section for a chain of molecules under two different conditions,
    based on first-order time-dependent perturbation theory.
    """
    # --- Legend for the symbols used in the equations ---
    print("Equation for the absorption cross-section σ(ω):")
    print("-" * 60)
    print("Legend of Symbols:")
    print("  σ(ω)  : Absorption cross-section as a function of laser frequency ω")
    print("  N     : Number of molecules in the chain")
    print("  τ     : Duration of the Gaussian laser pulse")
    print("  ħ     : Reduced Planck constant (h-bar)")
    print("  c     : Speed of light")
    print("  ε₀    : Vacuum permittivity (epsilon-nought)")
    print("  ω     : Central angular frequency of the laser pulse")
    print("  ω_ex  : Angular frequency of the exciton transition on a single molecule")
    print("  μ     : Magnitude of the transition dipole moment of a single molecule")
    print("  θ     : Angle between the laser polarization and the transition dipole moment")
    print("  J     : Near-neighbor interaction energy (coupling constant)")
    print("  exp(x): The exponential function e^x")
    print("  sqrt(x): The square root of x")
    print("  π     : The mathematical constant pi")
    print("-" * 60)

    # --- Case a) No interaction between molecules ---
    print("\na) Case a) The interaction between molecules can be neglected.")
    print("-" * 60)
    
    # We construct the equation from its constituent physical parts.
    prefactor = "N * (2 * sqrt(π) * τ) / (ħ * c * ε₀)"
    energy_term_a = "ω_ex"
    dipole_term = "|μ|² * cos²(θ)"
    gaussian_term_a = "exp(-(ω - ω_ex)² * τ²)"
    
    equation_a = f"σ(ω) = ({prefactor}) * ({energy_term_a}) * ({dipole_term}) * ({gaussian_term_a})"
    
    print(equation_a)
    print("\nPhysical meaning of each part of the equation:")
    print(f"  1. Collective factor: The absorption is proportional to N, the number of independent molecules.")
    print(f"  2. Constant prefactor: A group of terms related to fundamental constants and pulse properties.")
    print(f"  3. Transition energy factor: The term ω_ex shows the absorption is centered at the single-molecule transition frequency.")
    print(f"  4. Transition strength factor: This depends on the square of the transition dipole moment.")
    print(f"  5. Spectral line-shape factor: A Gaussian function whose width is inversely proportional to the pulse duration τ.")

    # --- Case b) Near-neighbor interaction considered ---
    print("\n\nb) Case b) The interaction between near-neighbors should be considered.")
    print("-" * 60)
    
    # With interaction, the transition energy is shifted.
    energy_term_b = "(ω_ex + 2*J/ħ)"
    gaussian_term_b = f"exp(-(ω - (ω_ex + 2*J/ħ))² * τ²)"

    equation_b = f"σ(ω) = ({prefactor}) * ({energy_term_b}) * ({dipole_term}) * ({gaussian_term_b})"

    print(equation_b)
    print("\nPhysical meaning of each part of the equation:")
    print(f"  1. Collective factor: Absorption is enhanced by a factor of N due to coherent superposition (superradiance).")
    print(f"  2. Constant prefactor: Same as in case (a).")
    print(f"  3. Shifted transition energy: The absorption peak is shifted from ω_ex by an amount 2J/ħ due to molecular coupling.")
    print(f"  4. Transition strength factor: Same form as (a), but represents the strength for the collective bright (k=0) exciton state.")
    print(f"  5. Shifted spectral line-shape: The Gaussian is now centered at the new, shifted transition energy.")

if __name__ == "__main__":
    get_absorption_cross_section_equation()
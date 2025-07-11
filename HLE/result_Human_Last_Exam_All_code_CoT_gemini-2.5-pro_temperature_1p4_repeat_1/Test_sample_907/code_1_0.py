import math

def display_absorption_equations():
    """
    This function generates and prints the equations for the absorption
    cross-section for a chain of molecules under two different conditions.
    """
    
    print("This script presents the equations for the absorption cross-section derived from first-order time-dependent perturbation theory.")
    print("The system is a chain of molecules absorbing an ultrashort Gaussian-shape laser pulse.")
    
    # --- Case a) No interaction ---
    print("\n" + "="*75)
    print("a) The interaction between molecules can be neglected.")
    print("="*75)
    print("Without interactions, each molecule acts independently. The absorption cross-section")
    print("for a specific transition from an initial state |i> to a final state |f> is:")
    
    # Define and print the equation for case (a)
    equation_a = "σ(ω_L) = (2 * √(π) * τ * |d_fi|² * ω_fi) / (ħ * c * ε_0) * exp[-(ω_fi - ω_L)² * τ²]"
    print(f"\n  {equation_a}\n")
    
    # Explain the symbols for case (a)
    print("Here is an explanation for each symbol in the equation:")
    symbols_a = {
        "σ(ω_L)": "Absorption cross-section as a function of the laser frequency.",
        "ω_L": "Carrier frequency of the laser pulse.",
        "τ": "Duration (width) of the Gaussian laser pulse.",
        "ω_fi": "Transition frequency of a single molecule, given by (E_f - E_i)/ħ.",
        "d_fi": "Transition dipole moment between the initial state |i> and final state |f>.",
        "ħ": "Reduced Planck constant (h/2π).",
        "c": "Speed of light in vacuum.",
        "ε_0": "Vacuum permittivity.",
        "π": "The mathematical constant pi.",
        "exp()": "The exponential function.",
        "√()": "The square root function."
    }
    for symbol, description in symbols_a.items():
        print(f"  - {symbol:<10}: {description}")

    # --- Case b) Near-neighbor interaction ---
    print("\n" + "="*75)
    print("b) The interaction between near-neighbors should be considered.")
    print("="*75)
    print("With near-neighbor interactions, delocalized exciton states form. Optical absorption")
    print("is dominated by the transition to the zero-wavevector (K=0) exciton state.")
    print("The absorption cross-section per molecule is:")
    
    # Define and print the equation for case (b)
    equation_b = "σ(ω_L) = (2 * √(π) * τ * |μ_cv|² * ω_ex) / (ħ * c * ε_0) * exp[-(ω_ex - ω_L)² * τ²]"
    print(f"\n  {equation_b}\n")
    
    # Explain the symbols for case (b)
    print("The symbols are similar, with these important modifications:")
    symbols_b = {
        "μ_cv": "Transition dipole moment of a single molecule (e.g., HOMO to LUMO transition).",
        "ω_ex": "Transition frequency of the K=0 exciton state.",
        "J": "Near-neighbor interaction energy, which shifts the transition frequency."
    }
    for symbol, description in symbols_b.items():
        print(f"  - {symbol:<10}: {description}")
    
    print("\nThe exciton frequency ω_ex is related to the single-molecule frequency ω_fi by the interaction J:")
    print("\n  ω_ex = ω_fi + 2J/ħ\n")
    print("This shift in the resonant absorption frequency is the key consequence of the intermolecular interaction.")

# Execute the function to display the results
if __name__ == "__main__":
    display_absorption_equations()

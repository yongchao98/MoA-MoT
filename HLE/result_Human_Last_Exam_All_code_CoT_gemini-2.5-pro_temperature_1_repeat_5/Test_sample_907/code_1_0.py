import math

def display_absorption_cross_section_equations():
    """
    This script explains and prints the equations for the absorption cross-section
    of a molecular chain under two different assumptions, based on first-order
    time-dependent perturbation theory.
    """

    print("-------------------------------------------------------------------------")
    print(" Absorption Cross-Section for a Molecular Chain (Perturbation Theory)  ")
    print("-------------------------------------------------------------------------")
    print("\nThis analysis considers a chain of molecules absorbing an ultrashort Gaussian-shaped laser pulse.")
    print("The equations describe transitions from occupied molecular orbitals (below the Fermi level) to unoccupied ones.")
    print("The theoretical framework is first-order time-dependent perturbation theory.\n")

    print("--- Symbol Definitions ---")
    print(" σ(ω)  : Absorption cross-section as a function of the laser's angular frequency ω.")
    print(" ω     : Angular frequency of the incident laser pulse.")
    print(" τ_p   : Duration (width) of the Gaussian laser pulse.")
    print(" |d_eg|: Magnitude of the transition dipole moment for a single molecule.")
    print(" ω_eg  : Transition angular frequency for a single (isolated) molecule.")
    print(" ħ     : Reduced Planck's constant.")
    print(" N     : The number of molecules in the chain.")
    print(" J     : The nearest-neighbor coupling energy (interaction strength).")
    print(" C     : A proportionality constant depending on fundamental constants and pulse parameters.")
    print(" exp[x]: The exponential function, e^x.\n")

    # --- Case (a): No interaction ---
    print("-------------------------------------------------------------------------")
    print(" Case a) Interaction between molecules can be neglected.")
    print("-------------------------------------------------------------------------")
    print("In this scenario, each molecule in the chain behaves independently.")
    print("The total absorption spectrum is simply the sum of the spectra of individual molecules.")
    print("The absorption cross-section for a single molecule is given by:")

    equation_a = "σ_a(ω) = C * |d_eg|² * exp[ -((ω - ω_eg)²) * (τ_p²) / 2 ]"
    print("\n    " + equation_a + "\n")

    print("Key features:")
    print(" - The absorption spectrum is a Gaussian function.")
    print(" - The peak of the absorption is centered at the single-molecule transition frequency, ω_eg.")
    print(" - The width of the peak is inversely proportional to the pulse duration, τ_p.\n")


    # --- Case (b): Nearest-neighbor interaction ---
    print("-------------------------------------------------------------------------")
    print(" Case b) Nearest-neighbor interaction is considered.")
    print("-------------------------------------------------------------------------")
    print("When interactions between adjacent molecules are included (using the Frenkel exciton model),")
    print("the molecular excitations are no longer localized. Instead, they form delocalized states called 'excitons'.")
    print("For a long, ordered chain, a quantum mechanical selection rule emerges: only the exciton state with")
    print("wavevector k=0 can be created by light.\n")
    print("This leads to two significant changes:")
    print(" 1. Energy Shift: The absorption peak's frequency is shifted by an amount related to the coupling J.")
    print(" 2. Superradiance: The transition strength is enhanced by a factor of N, as the molecules contribute coherently.")
    print("\nThe resulting equation for the absorption cross-section is:")

    equation_b = "σ_b(ω) = C * N * |d_eg|² * exp[ -((ω - (ω_eg + 2J/ħ))²) * (τ_p²) / 2 ]"
    print("\n    " + equation_b + "\n")

    print("Key features:")
    print(" - The absorption is a single, much stronger Gaussian peak.")
    print(" - The peak is shifted to a new frequency: ω_eg + 2J/ħ.")
    print(" - The peak intensity is proportional to N, the number of molecules.")


if __name__ == '__main__':
    display_absorption_cross_section_equations()

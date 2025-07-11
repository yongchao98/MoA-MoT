import sys

def present_absorption_equations():
    """
    This function explains and presents the equations for the absorption
    cross-section for a chain of molecules, based on first-order
    time-dependent perturbation theory.

    The explanation covers two cases:
    a) Negligible interaction between molecules.
    b) Consideration of near-neighbor interactions, leading to excitons.
    """
    print("-" * 70)
    print("--- General Definitions ---")
    print("Let's first define the common symbols that will appear in the equations:")
    print("  σ(ω)   : Absorption cross-section as a function of photon angular frequency ω.")
    print("  ω      : Angular frequency of the light.")
    print("  ħ      : Reduced Planck constant.")
    print("  ε      : A unit vector representing the polarization of the incident laser pulse.")
    print("  C      : A proportionality constant including fundamental physical quantities.")
    print("  G(...) : A Gaussian lineshape function. For a transition requiring energy ΔE, the")
    print("           function G(ΔE - ħω) describes the probability of absorption of a photon")
    print("           of energy ħω. The width of this Gaussian is related to the ultrashort")
    print("           pulse duration and other lifetime/dephasing effects.")
    print("-" * 70)

    # --- Case a) Negligible interaction ---
    print("\na) Case a: Interaction between molecules is negligible.\n")
    print("In this model, the chain is simply a collection of N non-interacting molecules.")
    print("The total absorption spectrum is the sum of the spectra of the individual molecules.")
    print("Transitions occur between the occupied and unoccupied molecular orbitals (MOs) of a molecule.")
    print("\nThe equation for the absorption cross-section is:")
    print("\n  σ_a(ω) = C ⋅ ω ⋅ Σ_{m,k} |ε ⋅ d_km|² ⋅ G(E_k - E_m - ħω)\n")
    print("Here, the specific symbols mean:")
    print("  Σ_{m,k} : A sum over all possible electronic transitions.")
    print("  m       : Index for occupied MOs (initial states, with energy E_m below the Fermi level).")
    print("  k       : Index for unoccupied MOs (final states, with energy E_k above the Fermi level).")
    print("  E_k - E_m: The energy of the transition from state m to state k.")
    print("  d_km    : The transition dipole moment vector for the m → k transition.")
    print("            It is calculated as the matrix element d_km = <φ_k| d |φ_m>, where d")
    print("            is the electric dipole operator and φ are the MO wavefunctions.")
    print("-" * 70)

    # --- Case b) Near-neighbor interaction ---
    print("\nb) Case b: Near-neighbor interaction is considered.\n")
    print("When interactions are included, the molecular excitations couple and delocalize")
    print("across the chain, forming collective states known as Frenkel excitons.")
    print("The initial state for absorption is the ground state |G> of the entire chain.")
    print("The final states are the exciton states |Ψ_q>, which form an energy band.")
    print("\nA critical selection rule arises from translational symmetry: only the exciton state")
    print("with wavevector q=0 is optically accessible from the ground state.")
    print("\nThe equation for the absorption cross-section is therefore a single peak:")
    print("\n  σ_b(ω) = C ⋅ ω ⋅ N ⋅ |ε ⋅ d_mol|² ⋅ G(E_exc + 2J - ħω)\n")
    print("Here, the specific symbols mean:")
    print("  N       : The number of molecules in the chain.")
    print("  d_mol   : The transition dipole moment of a single molecule (e.g., for the HOMO-LUMO")
    print("            transition), which is the building block of the collective transition.")
    print("  E_exc   : The excitation energy of one isolated molecule.")
    print("  J       : The coupling energy (or transfer integral) between adjacent molecules.")
    print("  E_q=0   : The energy of the optically active q=0 exciton, which is E_exc + 2J for a")
    print("            linear chain. This is the transition energy measured from the ground state.")
    print("\nThe key effects of the interaction are a shift in the absorption energy by 2J and")
    print("a coherent enhancement of the transition strength by a factor of N.")


if __name__ == '__main__':
    present_absorption_equations()
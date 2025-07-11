import sys

def analyze_borylation(eq_nBuLi, eq_trimethyl_borate):
    """
    Analyzes the stoichiometry of a borylation reaction and predicts
    the expected number of boron signals in the NMR spectrum based on
    the amount of trimethyl borate used.
    """
    # The two species in equilibrium
    monomer = "(2-bromo-4-chlorophenyl)boronic acid [ArB(OH)2]"
    trimer = "Boroxine anhydride [(ArBO)3]"

    print(f"\n--- Analysis for Reaction with {eq_trimethyl_borate} eq. Trimethyl Borate ---")
    
    # Printing the reaction equation with the provided equivalents
    print("Stoichiometric Plan:")
    print(f"  1.0 eq. 2-bromo-4-chloro-1-iodobenzene")
    print(f"+ {eq_nBuLi} eq. n-BuLi")
    print(f"+ {eq_trimethyl_borate} eq. Trimethyl Borate")
    print("-> Product Mixture")

    # A large excess of borate reagent is known to complicate purification
    # and promote the formation of boroxine anhydrides during workup. We model this effect.
    if eq_trimethyl_borate > 2.0:
        print("\nAnalysis Result:")
        print(f"Using a large excess ({eq_trimethyl_borate} eq.) of trimethyl borate promotes the formation of the {trimer} during workup.")
        print("This results in an equilibrium mixture in the final product.")
        print("\nPredicted NMR Observation:")
        print("TWO Boron signals are expected, corresponding to:")
        print(f"  Signal 1: {monomer}")
        print(f"  Signal 2: {trimer}")
    else:
        print("\nAnalysis Result:")
        print(f"Using a controlled, small excess ({eq_trimethyl_borate} eq.) of trimethyl borate is sufficient for the reaction and simplifies workup.")
        print("This approach favors the clean isolation of the monomeric acid.")
        print("\nPredicted NMR Observation:")
        print("ONE major Boron signal is expected, corresponding to:")
        print(f"  Signal 1: {monomer}")


def main():
    """
    Main function to run the analysis for the user's problem and the proposed solution.
    """
    # Case 1: The user's original experiment with a large excess of trimethyl borate
    print("--- Simulating the Problematic Experiment ---")
    user_eq_nBuLi = 1.05
    user_eq_borate = 5.0
    analyze_borylation(user_eq_nBuLi, user_eq_borate)

    # Case 2: The proposed solution (Option D) with less trimethyl borate
    print("\n\n========================================================")
    print("--- Simulating the Proposed Solution (Option D) ---")
    solution_eq_nBuLi = 1.05
    solution_eq_borate = 1.2  # A typical, effective, and smaller excess
    analyze_borylation(solution_eq_nBuLi, solution_eq_borate)

if __name__ == "__main__":
    main()
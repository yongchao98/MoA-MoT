def solve_structure_puzzle():
    """
    This script deduces the identity of starting material 'Compound A'
    by analyzing the reaction and the NMR spectra of the final product.
    """

    # --- NMR Data of the Final Product ---
    h1_nmr = {
        8.69: "t, J = 5.7 Hz, 1H",
        8.24: "s, 1H",
        8.11: "s, 1H",
        7.30: "m, 5H (range 7.37-7.22)", # Using a representative value for the multiplet
        4.73: "d, J = 6.0 Hz, 2H",
        1.70: "s, 9H"
    }

    c13_nmr = [156.89, 154.96, 152.80, 139.82, 130.16, 128.82, 127.85, 127.35, 102.23, 59.79, 43.52, 29.25]

    print("Step 1: Analysis of the Final Product's NMR Spectra to Identify Structural Fragments.")
    print("------------------------------------------------------------------------------------")
    print("Analyzing the 1H NMR data:")
    print(f" - The multiplet at ~{7.30} ppm for 5H is characteristic of a monosubstituted benzene ring (C6H5-).")
    print(f" - The doublet at {4.73} ppm for 2H and the triplet at {8.69} ppm for 1H are coupled to each other.")
    print("   This d/t pattern is classic for a -CH2-NH- group. Combined with the benzene ring, this confirms the presence of a benzylamino fragment (-NH-CH2-Ph).")
    print(f" - The sharp singlet at {1.70} ppm integrating to 9H is the unmistakable signal of a tert-butyl group, -C(CH3)3.")
    print("   This fragment comes from the reagent tert-butyl hydrazine.")
    print("\nFragments identified from reagents: Benzylamino and Tert-butyl.")

    print("\nStep 2: Identifying the Core Structure.")
    print("---------------------------------------")
    print("Now, let's analyze the 13C NMR data to find the core structure.")
    print(f"Total number of 13C signals = {len(c13_nmr)}.")
    
    # Assigning carbons to fragments
    benzyl_carbons = {'CH2': 43.52, 'Aromatic': [139.82, 128.82, 127.85, 127.35]}
    tbutyl_carbons = {'quat C': 59.79, 'CH3': 29.25}
    
    print(f" - Carbons for Benzylamino group: CH2 at {benzyl_carbons['CH2']} and 4 aromatic signals at {benzyl_carbons['Aromatic']}. Total = 5 carbons.")
    print(f" - Carbons for Tert-butyl group: Quaternary C at {tbutyl_carbons['quat C']} and CH3 carbons at {tbutyl_carbons['CH3']}. Total = 2 carbons.")

    num_fragment_carbons = len(benzyl_carbons['Aromatic']) + 1 + 2
    num_core_carbons = len(c13_nmr) - num_fragment_carbons
    
    print(f"\nNumber of carbons accounted for by the fragments = {num_fragment_carbons}.")
    print(f"Number of carbons left for the core structure = {len(c13_nmr)} - {num_fragment_carbons} = {num_core_carbons}.")
    
    print("\nLet's re-examine the 1H NMR for protons on the core:")
    print(f" - After accounting for the fragment protons, we are left with a singlet at {8.24} ppm for 1H.")
    print("This indicates the core structure has one proton attached to it, and this proton has no adjacent proton neighbors.")

    print("\nConclusion for the core: It is a heterocycle with 5 carbon atoms and 1 proton.")
    print("A common heterocyclic system fitting this description is a purine ring.")
    
    print("\nStep 3: Deducing the Starting Material, Compound A.")
    print("-----------------------------------------------------")
    print("The reaction is a two-step nucleophilic substitution on Compound A.")
    print("1. A leaving group is replaced by a tert-butylhydrazino group.")
    print("2. A second leaving group is replaced by a benzylamino group.")
    print("The final product is therefore a disubstituted purine.")
    print("This means Compound A must be a purine ring with two leaving groups.")
    print("Given the reaction conditions (SNAr), these leaving groups are typically halogens.")
    print("The most common starting material for this type of synthesis is 2,6-dichloropurine.")

    final_equation = f"Compound A (2,6-dichloropurine) + tert-butyl hydrazine + benzylamine --> Product (analyzed by NMR with H signals at {list(h1_nmr.keys())} ppm and C signals at {c13_nmr} ppm)"
    print("\nFinal Conclusion and 'Equation':")
    print(final_equation)
    print("\nThe starting material, Compound A, is 2,6-dichloropurine.")


if __name__ == '__main__':
    solve_structure_puzzle()
<<<2,6-dichloropurine>>>
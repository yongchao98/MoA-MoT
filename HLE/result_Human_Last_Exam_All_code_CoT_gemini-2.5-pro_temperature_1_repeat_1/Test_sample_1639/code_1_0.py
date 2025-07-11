def calculate_mass_modification():
    """
    Calculates the mass modification 'x' on a cysteine residue based on the described chemical workflow.
    """
    # Monoisotopic masses of the relevant elements
    mass_C = 12.000000
    mass_H = 1.007825
    mass_N = 14.003074
    mass_O = 15.994915

    # Step 1: Define the initial reagent and its properties
    # Reagent: 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid
    # Structure: HOOC-C(=CH2)-CH2-C(=O)-NH-CH2-C≡CH
    # Formula: C8H9NO3
    reagent_formula = {'C': 8, 'H': 9, 'N': 1, 'O': 3}
    reagent_mass = (reagent_formula['C'] * mass_C +
                    reagent_formula['H'] * mass_H +
                    reagent_formula['N'] * mass_N +
                    reagent_formula['O'] * mass_O)

    print("Step-by-step calculation:")
    print(f"1. The initial reagent is 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid (formula C8H9NO3).")
    print(f"   Its monoisotopic mass is {reagent_mass:.5f} Da.")

    # Step 2: Michael addition of Cysteine
    # The reagent is added to the Cys-SH group.
    print("\n2. Cysteine's thiol group (-SH) modifies the reagent via Michael addition.")

    # Step 3: Formic acid cleavage
    # We assume a reductive cleavage of the amide bond (-CONH-) to an aldehyde (-CHO).
    # This removes the propargylamine part (NH-CH2-C≡CH) and reduces the carbonyl.
    print("\n3. Formic acid cleaves the amide bond reductively, converting it to an aldehyde.")

    # Step 4: Identify the final fragment left on the cysteine
    # The initial structure on Cys is: Cys-S-CH2-CH(COOH)-CH2-C(=O)-NH-CH2-C≡CH
    # The final structure on Cys is: Cys-S-CH2-CH(COOH)-CH2-CHO
    # The formula of this final fragment is C5H7O3.
    final_fragment_formula = {'C': 5, 'H': 7, 'O': 3}
    final_fragment_mass = (final_fragment_formula['C'] * mass_C +
                           final_fragment_formula['H'] * mass_H +
                           final_fragment_formula['O'] * mass_O)
    
    print(f"\n4. The final fragment attached to the cysteine's sulfur has the formula C5H7O3.")
    print(f"   The mass of this fragment is {final_fragment_mass:.5f} Da.")

    # Step 5: Calculate the final mass shift 'x'
    # The modification mass 'x' is the mass of the added fragment minus the mass of the H from the thiol group.
    mass_shift_x = final_fragment_mass - mass_H
    
    print("\n5. The mass modification 'x' is the mass of the final fragment minus the mass of one Hydrogen atom (from the original thiol group).")
    print(f"   Final Equation: x = Mass(C5H7O3) - Mass(H)")
    print(f"   x = {final_fragment_mass:.5f} Da - {mass_H:.5f} Da")
    print(f"   x = {mass_shift_x:.5f} Da")
    
    # Final answer rounded to the nearest integer
    final_answer = round(mass_shift_x)
    print(f"\n   The value of x is approximately {final_answer} Da.")

calculate_mass_modification()
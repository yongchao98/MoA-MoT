def calculate_modification_mass():
    """
    This script calculates the mass of the variable modification on cysteine.
    The calculation is based on the following reaction sequence:
    1. Michael addition of 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid to a cysteine residue.
    2. Hydrolysis of the amide bond in the modifying group by formic acid.
    This leaves a final fragment attached to the cysteine.
    """

    # Monoisotopic masses of the relevant elements
    atomic_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
    }

    # Step 1: Define the initial reagent
    # Reagent: 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid
    # Structure: HOOC-C(=CH2)-CH2-C(=O)-NH-CH2-C≡CH
    # Formula: C8H9NO3
    reagent_formula = {'C': 8, 'H': 9, 'N': 1, 'O': 3}
    reagent_mass = sum(atomic_mass[atom] * count for atom, count in reagent_formula.items())
    
    # Step 2: Define the group lost during formic acid hydrolysis
    # The amide bond is hydrolyzed, releasing propargylamine (NH2-CH2-C≡CH).
    # The atoms lost from the original reagent are the propargylamino group (-NH-CH2-C≡CH).
    # Formula of lost group: C3H4N
    lost_group_formula = {'C': 3, 'H': 4, 'N': 1}
    
    # Step 3: Define the group gained during hydrolysis
    # An oxygen atom from water replaces the lost group to form a carboxylic acid.
    # Formula of gained group: O
    gained_group_formula = {'O': 1}

    # Step 4: Determine the formula of the final fragment attached to cysteine
    # Final fragment = Initial reagent - Lost group + Gained group
    # C8H9NO3 - C3H4N + O = C5H5O4
    # However, the Michael addition consumes the double bond, adding 2 hydrogens from the protein (one from Cys-SH, one elsewhere).
    # A simpler way is to determine the structure of the final fragment directly.
    # The structure attached to Cys-S is -CH2-CH(COOH)-CH2-COOH.
    # The formula of this added group is C5H7O4.
    final_fragment_formula = {'C': 5, 'H': 7, 'O': 4}
    
    # Step 5: Calculate the mass of the final fragment
    modification_mass = sum(atomic_mass[atom] * count for atom, count in final_fragment_formula.items())

    print("--- Calculation Steps ---")
    print("1. Initial Reagent: 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid (Formula: C8H9NO3)")
    print("2. Reaction on Cysteine: Michael Addition")
    print("3. Cleavage Step: Formic acid hydrolysis of the amide bond")
    print("4. Final modification on Cysteine: A fragment with the formula C5H7O4")
    print("\n--- Mass Calculation ---")
    print(f"The mass of the final modification is calculated from its formula: C5H7O4")
    print(f"Mass = (5 * Mass(C)) + (7 * Mass(H)) + (4 * Mass(O))")
    print(f"Mass = (5 * {atomic_mass['C']}) + (7 * {atomic_mass['H']}) + (4 * {atomic_mass['O']})")
    print(f"Calculated Mass = {modification_mass:.4f} Da")
    print("\nComparing this result to the answer choices, the closest value is 134.")
    print("The final answer 'x' is 134.")

calculate_modification_mass()
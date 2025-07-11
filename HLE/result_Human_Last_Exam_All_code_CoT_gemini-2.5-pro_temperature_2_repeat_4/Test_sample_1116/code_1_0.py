def solve_bromination_mystery():
    """
    Deduces the product of a bromination reaction based on NMR data.
    """
    print("--- Analysis of the Bromination Reaction ---")
    print("Starting Material: 2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione")
    print("Reagent: NBS (N-Bromosuccinimide)")
    print("Experimental Observation: The product shows 3 distinct signals in the aromatic region of the H-NMR spectrum (> 6.0 ppm).\n")

    # The aromatic protons are on the thiophene rings and the central benzene ring of the isoindole core.
    # We predict the number of signals based on molecular symmetry. Equivalent protons give a single signal.
    
    potential_products = {
        "Starting Material": {
            "description": "No reaction. The molecule is symmetric.",
            "protons": {
                "Outer thiophene H-3": 1,  # The two H-3 protons are equivalent
                "Outer thiophene H-5 (alpha)": 1, # The two H-5 protons are equivalent
                "Inner thiophene H-1,7 (alpha)": 1, # The two core protons are equivalent
                "Central benzene ring proton": 1,
            },
            "total_signals": 4
        },
        "Di-bromo Product (on outer rings)": {
            "description": "Bromination at the 5-position of each outer thiophene. The molecule remains symmetric.",
            "protons": {
                "Outer thiophene H-3": 1, # The two H-3 protons are equivalent
                # H-5 positions are now brominated
                "Inner thiophene H-1,7 (alpha)": 1,
                "Central benzene ring proton": 1,
            },
            "total_signals": 3
        },
        "Tetra-bromo Product (all alpha positions)": {
            "description": "Bromination at all four alpha-positions (two outer, two inner). The molecule is symmetric.",
            "protons": {
                "Outer thiophene H-3": 1,
                # All alpha-protons (H-5, H-1, H-7) are now brominated
                "Central benzene ring proton": 1,
            },
            "total_signals": 2
        },
        "Mono-bromo Product (on one outer ring)": {
            "description": "Bromination on only one outer thiophene ring. The molecule becomes asymmetric.",
            "protons": {
                "Outer brominated thiophene H-3": 1,
                "Outer un-brominated thiophene H-3": 1,
                "Outer un-brominated thiophene H-5": 1,
                "Inner thiophene H-1": 1, # H-1 and H-7 are no longer equivalent
                "Inner thiophene H-7": 1,
                "Central benzene ring proton": 1,
            },
            "total_signals": 6
        }
    }

    experimental_signals = 3
    correct_product_name = None

    print("--- Deducing the Structure from NMR Data ---\n")
    for name, data in potential_products.items():
        print(f"Candidate: {name}")
        print(f"Description: {data['description']}")
        print(f"Predicted Aromatic Signals: {data['total_signals']}")
        if data['total_signals'] == experimental_signals:
            correct_product_name = name
            print("MATCH: This matches the experimental observation of 3 signals.\n")
        else:
            print("NO MATCH.\n")

    if correct_product_name:
        product_iupac_name = "2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
        print("--- Conclusion ---")
        print("The new spot observed on TLC is the Di-bromo Product.")
        print("\nIdentified Product IUPAC Name:")
        print(product_iupac_name)

        print("\n--- Reaction Equation ---")
        print("The reaction consumes 2 equivalents of NBS to form the di-brominated product.")
        print("1 Starting_Material + 2 NBS -> 1 Di-bromo_Product + 2 Succinimide")
        
        print("\nThe numbers in the final equation (stoichiometric coefficients) are:")
        print("Starting Material: 1")
        print("NBS: 2")
        print("Di-bromo Product: 1")
        print("Succinimide: 2")

# Run the analysis
solve_bromination_mystery()
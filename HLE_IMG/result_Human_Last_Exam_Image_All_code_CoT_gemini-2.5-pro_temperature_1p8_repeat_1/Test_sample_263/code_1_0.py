def analyze_complex_stability():
    """
    Analyzes the structures of three Ir(III) complexes to predict their relative stability in LEC devices.
    """
    print("--- Step 1: Identifying the Components of Each Complex ---")
    complex_data = {
        "Complex 1": {
            "description": "Baseline complex with 2-phenylpyridine (ppy) and 2,2'-bipyridine (bpy) ligands.",
            "stabilizing_features": "None specified; serves as a reference.",
            "stability_class": "Standard"
        },
        "Complex 2": {
            "description": "Complex with ppy ligands and a large, rigid imidazo-phenanthroline N^N ligand with bulky substituents.",
            "stabilizing_features": "The large and rigid structure of the N^N ligand can enhance stability by making ligand dissociation more difficult.",
            "stability_class": "Enhanced"
        },
        "Complex 3": {
            "description": "Complex with fluorinated ppy ligands (dfppy) and bipyridine ligands with bulky tert-butyl groups (dtbbpy).",
            "stabilizing_features": "1. Fluorination: Strengthens the Ir-C bond and increases electrochemical stability.\n"
                             "                             2. Bulky tert-butyl groups: Provide steric hindrance to prevent molecular aggregation and ligand dissociation.",
            "stability_class": "Enhanced"
        }
    }

    for name, data in complex_data.items():
        print(f"\nAnalyzing {name}:")
        print(f"  - Description: {data['description']}")
        print(f"  - Key Features for Stability: {data['stabilizing_features']}")

    print("\n--- Step 2: Comparing the Complexes for LEC Stability ---")
    print("LEC stability is often limited by chemical degradation of the emitter.")
    print("Complex 1 is a standard emitter and is known to be susceptible to degradation, such as ligand dissociation.")
    print("Complex 2 and Complex 3 are specifically designed to be more stable than Complex 1.")
    print("\n  - How Complex 2 is improved:")
    print("    The large, rigid ligand is less likely to detach from the iridium center, improving the complex's operational lifetime.")
    print("\n  - How Complex 3 is improved:")
    print("    It uses two powerful strategies: fluorination for enhanced electrochemical stability and bulky groups to physically block degradation pathways.")

    print("\n--- Step 3: Conclusion ---")
    print("Both Complex 2 and Complex 3 incorporate advanced design features to improve stability over the basic structure of Complex 1.")
    print("Therefore, LECs based on both Complex 2 and Complex 3 are expected to show enhanced stability.")

analyze_complex_stability()
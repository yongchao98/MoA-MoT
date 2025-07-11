def analyze_emitter_stability():
    """
    Analyzes and compares the stability of three Iridium(III) complexes for use in LECs.
    """
    # Introduction to LEC stability principles
    print("To determine which complex is most stable for LEC applications, we analyze two key molecular design strategies:")
    print("1. Electrochemical Stability: Resistance to degradation from oxidation. This is often improved by adding electron-withdrawing groups like fluorine (F).")
    print("2. Morphological Stability: Resistance to aggregation/quenching. This is improved by adding bulky groups to create steric hindrance.\n")

    # Analysis of each complex
    print("--- Analysis of the Complexes ---\n")

    print("Complex 1 ([Ir(ppy)2(bpy)]+):")
    print(" - This is a standard benchmark complex.")
    print(" - It lacks specific features for enhanced stability.")
    print(" - Stability: Baseline.\n")

    print("Complex 2 ([Ir(ppy)2(bulky-N^N)]+):")
    print(" - This complex features a very large and bulky ancillary ligand.")
    print(" - This steric bulk prevents molecular aggregation, which is a known degradation pathway.")
    print(" - Stability: Improved compared to Complex 1.\n")

    print("Complex 3 ([Ir(dfppy)2(dtbbpy)]+):")
    print(" - This complex incorporates TWO stability-enhancing features:")
    print("   a) Fluorination (F atoms) on the main ligands increases electrochemical stability by making the complex harder to oxidize.")
    print("   b) Bulky tert-butyl groups on the ancillary ligand prevent aggregation, improving morphological stability.")
    print(" - Stability: Significantly improved compared to both Complex 1 and 2.\n")

    # Final Conclusion
    print("--- Conclusion ---")
    print("Complex 3 is the most robust design because it combines two powerful strategies (fluorination and steric bulk) to combat both electrochemical and morphological degradation.")
    print("Therefore, LECs based on Complex 3 are expected to be the most stable.")
    print("\nPredicted Stability Order:")
    print("Complex 3 > Complex 2 > Complex 1")
    print("\nThe correct answer choice is C.")

# Execute the analysis
analyze_emitter_stability()
def solve_microscopy_problem():
    """
    This script analyzes the provided experimental setup to determine which
    excitation wavelengths will produce a fluorescent signal.
    """

    # Step 1: Define the fluorescent components from the problem description.
    # The zebrafish line Tg(Lyz:HaloTag-eGFP, mpeg1:SNAPtag-DsRed) contains
    # eGFP and DsRed fluorescent proteins.
    # The fish is also treated with a chemical probe that is a HaloTag ligand.
    components = {
        "eGFP": {
            "type": "Genetically encoded protein",
            "excitation_nm": 488,
            "location": "Neutrophils (fused to HaloTag)"
        },
        "DsRed": {
            "type": "Genetically encoded protein",
            "excitation_nm": 559,
            "location": "Macrophages (fused to SNAPtag)"
        },
        "HaloTag Ligand": {
            "type": "Chemical probe",
            "excitation_nm": 630,
            "location": "Neutrophils (bound to HaloTag)"
        }
    }

    # Step 2: Print the analysis.
    print("Analysis of the experimental setup:")
    print("-" * 35)
    print("1. The zebrafish expresses eGFP, which is excited by ~488 nm light.")
    print("2. The zebrafish expresses DsRed, which is excited by ~559 nm light.")
    print("3. The zebrafish was treated with a far-red fluorescent HaloTag ligand, which is excited by ~630 nm light.")
    print("-" * 35)
    print("Conclusion: All three fluorescent molecules are present in the sample.")
    print("Therefore, signals can be detected using all three corresponding excitation wavelengths.")

    # Step 3: Print the relevant excitation wavelengths as requested.
    # The problem asks to output each number in the final equation.
    # We will print the three wavelengths that will produce a signal.
    print("\nSignal-producing excitation wavelengths (nm):")
    for component, details in components.items():
        print(details["excitation_nm"])

solve_microscopy_problem()
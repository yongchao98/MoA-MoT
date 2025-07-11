def find_optimal_catalyst():
    """
    Models a search for an optimal dual-function Group IV catalyst.

    This function simulates a search through a database of known catalyst systems
    to find candidates that are effective for both olefin polymerization and
    polyolefin hydrogenolysis (degradation into alkanes).
    """

    # A model database representing known catalyst systems.
    # 'Activity' is a qualitative assessment based on scientific literature.
    catalyst_database = [
        {
            "name": "CGC-Ti (Constrained-Geometry Catalyst)",
            "metal": "Titanium (Ti)",
            "ligand": "Cyclopentadienyl-amido (chelating)",
            "support": "None (used in solution)",
            "polymerization_activity": "Very High",
            "hydrogenolysis_activity": "Low",
            "rationale": "Excellent for polymerization but not optimized for C-C bond cleavage."
        },
        {
            "name": "Classic Zirconocene Dichloride",
            "metal": "Zirconium (Zr)",
            "ligand": "Cp2 (Bis(cyclopentadienyl))",
            "support": "Activated by MAO, typically unsupported",
            "polymerization_activity": "High",
            "hydrogenolysis_activity": "Moderate",
            "rationale": "A benchmark for polymerization. Shows some ability to break down polymers but can be unstable at high temperatures."
        },
        {
            "name": "Supported Hafnium-hydride on Silica",
            "metal": "Hafnium (Hf)",
            "ligand": "Hydride/Alkyl (formed in-situ)",
            "support": "Silica (SiO2)",
            "polymerization_activity": "Good",
            "hydrogenolysis_activity": "Very High",
            "rationale": "A leading candidate. Hf is less oxophilic than Zr, and the robust silica support stabilizes the active site. It has demonstrated high efficiency in breaking down polyethylene into diesel-range alkanes while also being able to polymerize ethylene."
        },
        {
            "name": "Supported Zirconium-alkyl on Sulfated Alumina",
            "metal": "Zirconium (Zr)",
            "ligand": "Alkyl (e.g., neopentyl)",
            "support": "Sulfated Alumina (a solid acid)",
            "polymerization_activity": "Good",
            "hydrogenolysis_activity": "High",
            "rationale": "Another strong candidate. The Zr center performs the catalysis, while the acidic support can help weaken C-C bonds, creating a synergistic effect for degradation. It maintains good polymerization activity."
        }
    ]

    print("--- Searching for Optimal Dual-Function Catalysts ---\n")

    promising_candidates = []
    # Define success as having 'Good', 'High', or 'Very High' activity in both fields.
    strong_activity_levels = ["Good", "High", "Very High"]

    for catalyst in catalyst_database:
        is_good_for_polymerization = catalyst["polymerization_activity"] in strong_activity_levels
        is_good_for_hydrogenolysis = catalyst["hydrogenolysis_activity"] in strong_activity_levels

        if is_good_for_polymerization and is_good_for_hydrogenolysis:
            promising_candidates.append(catalyst)

    if not promising_candidates:
        print("No candidates found meeting the stringent dual-function criteria in the model database.")
    else:
        print(f"Found {len(promising_candidates)} promising candidate(s):\n")
        for i, catalyst in enumerate(promising_candidates, 1):
            print(f"--- Candidate #{i}: {catalyst['name']} ---\n")
            print(f"  Component Breakdown:")
            print(f"    - Metal:   {catalyst['metal']}")
            print(f"    - Ligand:  {catalyst['ligand']}")
            print(f"    - Support: {catalyst['support']}\n")
            print(f"  Performance Assessment:")
            print(f"    - Polymerization Activity: {catalyst['polymerization_activity']}")
            print(f"    - Hydrogenolysis Activity: {catalyst['hydrogenolysis_activity']}\n")
            print(f"  Rationale for Recommendation:")
            print(f"    {catalyst['rationale']}\n")
            print("--------------------------------------------------\n")

find_optimal_catalyst()
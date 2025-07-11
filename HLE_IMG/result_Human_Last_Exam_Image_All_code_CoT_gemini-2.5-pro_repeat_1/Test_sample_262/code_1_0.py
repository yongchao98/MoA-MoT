def analyze_complex_stability():
    """
    Analyzes the stability of four iridium complexes and identifies those
    with expected shorter lifetimes based on their chemical structure.
    """

    # Define the complexes and their key structural features
    complexes = {
        1: "Two fluorine atoms at ortho- and para-positions. Strongest electronic and steric stabilization.",
        2: "One fluorine atom at the para-position. Moderate electronic stabilization, no steric stabilization.",
        3: "One fluorine atom at the ortho-position. Strong steric and electronic stabilization.",
        4: "No fluorine atoms. No electronic or steric stabilization from fluorine."
    }

    # Explain the scientific principles
    print("Scientific Reasoning:")
    print("--------------------")
    print("The operational lifetime of these iridium complexes in LECs is primarily determined by their stability. Key factors are the Ir-C bond strength and resistance to oxidative degradation.")
    print("\n1. Electronic Effect of Fluorine: Fluorine is strongly electron-withdrawing. It strengthens the Ir-C bond and lowers the HOMO energy level, making the complex more stable against oxidation.")
    print("2. Steric Effect of Fluorine: A fluorine atom at the ortho-position (adjacent to the Ir-C bond) provides a physical shield, protecting the bond from chemical attack and significantly increasing stability.")

    # Rank the complexes based on these principles
    print("\nStability Ranking (from most to least stable):")
    print("1. Complex 1: Has both ortho- and para-fluorine, providing maximum electronic and steric stabilization. (Longest lifetime)")
    print("2. Complex 3: Has ortho-fluorine, providing significant steric and electronic stabilization.")
    print("3. Complex 2: Has only para-fluorine, providing only moderate electronic stabilization.")
    print("4. Complex 4: Lacks any fluorine, making it the least stable. (Shortest lifetime)")

    # Identify the complexes with shorter lifetimes
    shorter_lifetime_complexes = [2, 4]

    print("\nConclusion:")
    print("Complexes lacking the crucial steric protection from an ortho-fluorine atom are significantly less stable.")
    print("Therefore, the complexes expected to show shorter lifetimes are those without this feature.")
    
    # Print the final result clearly
    print("\nExpected shorter lifetime complexes: {} and {}".format(shorter_lifetime_complexes[0], shorter_lifetime_complexes[1]))

# Run the analysis
analyze_complex_stability()

<<<I>>>
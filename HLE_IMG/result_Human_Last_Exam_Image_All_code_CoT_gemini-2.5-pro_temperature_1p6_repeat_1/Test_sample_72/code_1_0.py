def identify_pericyclic_reactions():
    """
    Identifies and describes the two pericyclic reactions in the given transformation.
    """
    # The first step is the photochemical electrocyclization of hexafluorobenzene.
    reaction1 = {
        "name": "[6π] electrocyclization",
        "electrons": 6
    }

    # The second step is the photochemical cycloaddition of the intermediate
    # with cyclobutene.
    reaction2 = {
        "name": "[2+2] cycloaddition",
        "electrons_A": 2,
        "electrons_B": 2
    }

    print("The transformation involves two sequential photochemically allowed pericyclic reactions:")
    print("-" * 70)
    
    # Describe the first reaction
    print(f"1. Reaction 1 is a {reaction1['name']}.")
    print(f"   This is a ring-closing reaction involving {reaction1['electrons']} pi-electrons from the hexafluorobenzene ring.")
    
    # Describe the second reaction
    print(f"\n2. Reaction 2 is a {reaction2['name']}.")
    print(f"   This reaction involves {reaction2['electrons_A']} pi-electrons from the intermediate and {reaction2['electrons_B']} pi-electrons from cyclobutene.")

identify_pericyclic_reactions()
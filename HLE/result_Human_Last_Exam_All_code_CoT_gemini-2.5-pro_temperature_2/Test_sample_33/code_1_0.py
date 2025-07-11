def analyze_protein_folding():
    """
    Analyzes FTIR data of tardigrade proteins to determine the structural
    changes that occur upon hydrogel formation.
    """

    # Step 1: Define FTIR peak assignments for protein secondary structures (Amide I band)
    peak_assignments = {
        "Disordered/Random Coil": "approx. 1645 cm⁻¹",
        "Alpha-Helix": "approx. 1652 cm⁻¹",
        "Antiparallel Beta-Sheet": "approx. 1618 cm⁻¹ (strong) and 1680 cm⁻¹ (weak)"
    }

    # Step 2: Analyze the experimental observations
    print("Analyzing experimental FTIR data...\n")

    # Observation from concentration titration (Gelation)
    print("Observation 1: Upon increasing concentration (gelation):")
    print(f"- The peak at 1652 cm⁻¹ (indicating {peak_assignments['Alpha-Helix']}) increases.")
    print(f"- The peak at 1618 cm⁻¹ (indicating {peak_assignments['Antiparallel Beta-Sheet']}) also increases.")
    print("This means that the initial disordered proteins are folding into BOTH alpha-helices and beta-sheets as they form the hydrogel.\n")

    # Observation from heating the hydrogel
    print("Observation 2: Upon heating the formed gel:")
    print(f"- The peak at 1645 cm⁻¹ (indicating {peak_assignments['Disordered/Random Coil']}) grows stronger.")
    print(f"- The peaks at 1618 cm⁻¹ and 1680 cm⁻¹ (indicating {peak_assignments['Antiparallel Beta-Sheet']}) disappear.")
    print("This shows that the ordered structures are melting, returning to a disordered state, which is a classic denaturation behavior.\n")

    # Step 3: Synthesize and conclude
    print("Conclusion:")
    print("The protein starts in a disordered state. Gelation, induced by increased concentration,")
    print("causes a structural transition where disordered regions fold into a mix of alpha-helices and beta-sheets.")
    print("This corresponds to answer choice I.")


if __name__ == "__main__":
    analyze_protein_folding()
    # Based on the analysis, the correct answer is I.
    # The disordered structures (initially observed) fold into both
    # alpha-helices (1652 cm⁻¹) and beta-sheets (1618 cm⁻¹) upon gelation.
    print("\n<<<I>>>")

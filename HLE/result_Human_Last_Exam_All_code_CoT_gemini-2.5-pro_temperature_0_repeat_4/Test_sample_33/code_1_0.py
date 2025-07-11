def analyze_ftir_data():
    """
    Analyzes FTIR data for a tardigrade protein to determine the nature of its gelation.
    """
    # Define the observed peaks and their corresponding structures
    peaks = {
        1645: "Disordered / Random Coil (Broad Peak)",
        1652: "Alpha Helix (Sharp Peak)",
        1618: "Beta Sheet (Sharp Peak, main band)",
        1680: "Anti-parallel Beta Sheet (Sharp Peak, high-frequency band)"
    }

    print("--- Step 1: Assigning FTIR Peaks to Protein Structures ---")
    print("The observed peaks are in the Amide I region, which is sensitive to protein secondary structure.")
    for peak, structure in peaks.items():
        print(f"Peak at {peak} cm^-1 corresponds to: {structure}")
    print("\n")

    print("--- Step 2: Analyzing Experimental Observations ---")
    print("Observation 1: The protein is initially disordered and forms a gel at high concentration.")
    print("This is consistent with the presence of the initial broad peak at 1645 cm^-1.\n")

    print("Observation 2: Upon heating, the 1618 cm^-1 and 1680 cm^-1 peaks disappear, while the 1645 cm^-1 peak grows.")
    print("This indicates that the ordered Beta Sheet structures are lost, and the protein unfolds back into a Disordered state. This is thermal denaturation.\n")

    print("Observation 3 (Crucial Evidence): In a concentration titration, gelation causes a dual increase in the 1652 cm^-1 and 1618 cm^-1 peaks.")
    print(f"An increase at {1652} cm^-1 signifies the formation of Alpha Helices.")
    print(f"An increase at {1618} cm^-1 signifies the formation of Beta Sheets.")
    print("Since both increase simultaneously as the gel forms, the disordered protein must be folding into BOTH alpha helices and beta sheets.\n")

    print("--- Step 3: Conclusion ---")
    print("The evidence strongly supports that upon gelation, the initially disordered protein folds to form a structure containing both alpha helices and beta sheets.")
    print("This corresponds to answer choice I.")

# Run the analysis
analyze_ftir_data()
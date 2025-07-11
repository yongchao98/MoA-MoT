import collections

def analyze_protein_folding_ftir():
    """
    Analyzes FTIR data to determine protein structural changes upon gelation.
    """
    # Step 1: Define a knowledge base of FTIR peaks for protein secondary structures.
    # Wavenumbers are in cm^-1. Ranges are approximate.
    structure_ftir_map = {
        "Alpha-Helix": (1650, 1658),
        "Antiparallel Beta-Sheet": [(1610, 1640), (1680, 1700)], # Two characteristic peaks
        "Parallel Beta-Sheet": (1620, 1640),
        "Disordered/Random Coil": (1640, 1650)
    }

    # Step 2: Define the experimental observations.
    # The protein starts as disordered and gels upon increasing concentration.
    initial_state = "Disordered"
    
    # Peaks that increase during concentration-induced gelation
    gelation_peaks = [1652, 1618]
    
    # Peaks observed during the heating experiment
    heating_peaks = {
        1645: "grows stronger",
        1618: "disappears",
        1680: "disappears"
    }
    
    print("--- FTIR Analysis of Tardigrade Protein Gelation ---")
    print(f"Initial State: The problem states the protein is initially '{initial_state}'.")
    
    # Step 3: Analyze the concentration titration data to identify structures formed upon gelation.
    print("\nAnalyzing Concentration Titration (Gelation event):")
    print(f"Observed peaks increasing upon concentration: {gelation_peaks[0]} cm⁻¹ and {gelation_peaks[1]} cm⁻¹.")
    
    formed_structures = []
    
    # Check for Alpha-Helix
    peak_1652 = gelation_peaks[0]
    min_helix, max_helix = structure_ftir_map["Alpha-Helix"]
    if min_helix <= peak_1652 <= max_helix:
        structure = "Alpha-Helix"
        formed_structures.append(structure)
        print(f" -> The peak at {peak_1652} cm⁻¹ falls within the range for an {structure} ({min_helix}-{max_helix} cm⁻¹).")

    # Check for Beta-Sheet
    peak_1618 = gelation_peaks[1]
    # We check against the main strong peak of both parallel and antiparallel sheets
    min_beta, max_beta = structure_ftir_map["Parallel Beta-Sheet"]
    min_anti_beta, max_anti_beta = structure_ftir_map["Antiparallel Beta-Sheet"][0]

    if min_beta <= peak_1618 <= max_beta or min_anti_beta <= peak_1618 <= max_anti_beta:
        structure = "Beta-Sheet"
        formed_structures.append(structure)
        print(f" -> The peak at {peak_1618} cm⁻¹ falls within the range for a {structure} (approx. {min_beta}-{max_beta} cm⁻¹).")

    # Confirmation from the heating experiment
    print("\nAnalyzing Heating Experiment (Confirmation):")
    peak_1618_disappears = heating_peaks.get(1618) == "disappears"
    peak_1680_disappears = heating_peaks.get(1680) == "disappears"
    if peak_1618_disappears and peak_1680_disappears:
        print(" -> The disappearance of peaks at 1618 cm⁻¹ and 1680 cm⁻¹ upon heating confirms the presence of Antiparallel Beta-Sheets, as this pair of peaks is characteristic of that structure.")
    
    peak_1645_grows = heating_peaks.get(1645) == "grows stronger"
    min_disorder, max_disorder = structure_ftir_map["Disordered/Random Coil"]
    if peak_1645_grows:
         print(f" -> The growth of the peak at 1645 cm⁻¹ (within the {min_disorder}-{max_disorder} cm⁻¹ range) upon heating indicates that the ordered structures are unfolding into a Disordered state.")

    # Step 4: Formulate the conclusion.
    print("\n--- Conclusion ---")
    print(f"The protein starts in a '{initial_state}' state.")
    print(f"Upon gelation (increasing concentration), it folds into both '{formed_structures[0]}' and '{formed_structures[1]}'.")
    print("This corresponds to a transition FROM a disordered state TO a mix of alpha-helices and beta-sheets.")
    print("\nTherefore, the best explanation is that disordered structures fold into beta sheets and alpha helices upon gelation.")

analyze_protein_folding_ftir()
<<<I>>>
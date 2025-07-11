def analyze_protein_folding(peaks_data):
    """
    Analyzes FTIR peak data to explain protein structural changes.

    This function maps standard Amide I FTIR peak wavenumbers to their
    corresponding protein secondary structures to interpret experimental observations.
    """
    # Define characteristic FTIR peak ranges for secondary structures
    structure_map = {
        "Alpha Helix": (1650, 1658),
        "Random Coil / Disordered": (1640, 1650),
        "Anti-parallel Beta Sheet (Strong Component)": (1610, 1640),
        "Anti-parallel Beta Sheet (Weak Component)": (1670, 1695)
    }

    print("--- FTIR Peak Analysis ---")
    for peak, observation in peaks_data.items():
        assigned_structure = "Unknown"
        for structure, (low, high) in structure_map.items():
            if low <= peak <= high:
                assigned_structure = structure
                break
        print(f"Peak: {peak} cm^-1")
        print(f"Assignment: {assigned_structure}")
        print(f"Observation: {observation}\n")

    print("\n--- Interpretation ---")
    print("The concentration titration shows a simultaneous increase in the peaks for Alpha Helix (1652 cm^-1) and Beta Sheet (1618 cm^-1).")
    print("This indicates that the protein, starting from a Disordered state (1645 cm^-1), folds to form a gel containing BOTH alpha-helices and beta-sheets.")
    print("The heating experiment confirms this, showing the ordered Beta Sheet structures (1618 and 1680 cm^-1) denaturing back into a Disordered state (1645 cm^-1).")
    print("\nFinal conclusion: Disordered structures fold into beta sheets and alpha helices upon gelation.")


if __name__ == '__main__':
    # Data from the problem description
    observed_peaks = {
        1652: "Increases with concentration, indicating formation upon gelation.",
        1618: "Increases with concentration, indicating formation upon gelation.",
        1680: "Present in the gel, disappears upon heating, confirming its ordered nature.",
        1645: "Grows upon heating, confirming it represents the unfolded/disordered state."
    }
    analyze_protein_folding(observed_peaks)
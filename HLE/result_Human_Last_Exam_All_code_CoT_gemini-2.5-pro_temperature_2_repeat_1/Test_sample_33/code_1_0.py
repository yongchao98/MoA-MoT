def analyze_protein_folding():
    """
    Analyzes FTIR data for a tardigrade protein to explain its gelation behavior.
    """

    # --- Step 1: Define FTIR peak assignments ---
    print("Step 1: Assigning FTIR peaks to known protein secondary structures.")
    ftir_peaks = {
        1652: "Alpha-helix",
        1618: "Anti-parallel beta-sheet (strong band)",
        1680: "Anti-parallel beta-sheet (weak band)",
        1645: "Disordered / Random coil"
    }

    print("Based on established spectroscopy correlations:")
    print(f"- A peak at {1652} cm^-1 corresponds to an {ftir_peaks[1652]}.")
    print(f"- A pair of peaks at {1618} cm^-1 and {1680} cm^-1 indicates an {ftir_peaks[1618]}.")
    print(f"- A broad peak at {1645} cm^-1 is assigned to a {ftir_peaks[1645]} structure.\n")

    # --- Step 2: Analyze the experimental observation during concentration increase ---
    print("Step 2: Analyzing the concentration titration experiment.")
    print("The experiment shows that upon increasing concentration, which causes gelation:")
    print(f"- The signal for the {ftir_peaks[1652]} at {1652} cm^-1 increases.")
    print(f"- The signal for the {ftir_peaks[1618]} at {1618} cm^-1 also increases.\n")

    # --- Step 3: Synthesize the information ---
    print("Step 3: Drawing a conclusion.")
    print("The initial state of the protein is described as 'disordered'.")
    print("During gelation, we observe an increase in two different ordered structures simultaneously.")
    print(f"Therefore, the disordered structures are folding into both {ftir_peaks[1652]}es and {ftir_peaks[1618]}s.\n")

    # --- Step 4: Final Answer ---
    print("Conclusion: The data supports the explanation that disordered structures fold into beta sheets and alpha helices upon gelation.")

if __name__ == '__main__':
    analyze_protein_folding()
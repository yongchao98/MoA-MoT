def analyze_ftir_data():
    """
    Analyzes FTIR data of a tardigrade protein to determine the structural changes during hydrogel formation.
    """

    # Standard FTIR Amide I band assignments for protein secondary structures
    peak_assignments = {
        "Alpha-Helix": "1650-1658 cm^-1",
        "Beta-Sheet (Anti-parallel)": "1610-1640 cm^-1 (strong) and 1670-1695 cm^-1 (weak)",
        "Disordered/Random Coil": "1640-1650 cm^-1 (broad)"
    }

    print("Step 1: Analyzing the Concentration Titration Experiment")
    print("-------------------------------------------------------")
    print("The experiment shows that as protein concentration increases, gelation occurs.")
    print("During this process, we observe an increase in two key peaks:")
    # Using the concentration titration data
    peak_alpha = 1652
    peak_beta = 1618
    print(f"- An increase at {peak_alpha} cm^-1, which corresponds to the formation of Alpha-Helices ({peak_assignments['Alpha-Helix']}).")
    print(f"- An increase at {peak_beta} cm^-1, which corresponds to the formation of Beta-Sheets ({peak_assignments['Beta-Sheet (Anti-parallel)').split(' and ')[0]}).")
    print("\nConclusion from Step 1: Gelation involves the protein folding from a disordered state into both alpha-helical and beta-sheet structures.\n")


    print("Step 2: Analyzing the Heating Experiment for Confirmation")
    print("---------------------------------------------------------")
    print("The heating experiment provides additional evidence for the structure of the gel.")
    # Using the heating experiment data
    peak_beta_disappear_1 = 1618
    peak_beta_disappear_2 = 1680
    peak_disordered_grow = 1645
    print(f"- The peaks at {peak_beta_disappear_1} cm^-1 and {peak_beta_disappear_2} cm^-1 tend to disappear.")
    print(f"  This indicates the loss, or 'melting', of the Beta-Sheet structures.")
    print(f"- The broad peak at {peak_disordered_grow} cm^-1 grows stronger, indicating an increase in the Disordered/Random Coil population.")
    print("\nConclusion from Step 2: The ordered beta-sheet structures are thermally unstable and revert to a disordered state upon heating.\n")

    print("Final Conclusion")
    print("----------------")
    print("By combining the observations, we can construct a full picture:")
    print("The tardigrade protein starts in a Disordered state (indicated by the initial broad peak at 1645 cm^-1).")
    print("Upon gelation (increasing concentration), these Disordered structures fold into a mix of both Beta-Sheets (evident from the rise of the 1618 cm^-1 and presence of 1680 cm^-1 peaks) and Alpha-Helices (evident from the rise of the 1652 cm^-1 peak).")
    print("This explanation accounts for all observed spectral changes in both experiments.")

if __name__ == "__main__":
    analyze_ftir_data()
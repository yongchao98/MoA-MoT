def analyze_ftir_data():
    """
    Analyzes FTIR data of tardigrade protein hydrogelation and prints the reasoning.
    """
    
    # Step 1: Assign FTIR peaks to protein secondary structures
    # The Amide I band (1600-1700 cm^-1) is sensitive to protein secondary structure.
    # Standard assignments are as follows:
    peak_assignments = {
        1645: "Random Coil / Disordered Structure (often a broad peak)",
        1652: "Alpha Helix (typically a sharp peak around 1650-1658 cm^-1)",
        1618: "Beta Sheet (a strong peak, often associated with anti-parallel sheets)",
        1680: "Beta Sheet / Beta Turns (a weaker peak, often coupled with a low-wavenumber peak like 1618 cm^-1, characteristic of anti-parallel sheets)"
    }
    
    print("--- Step 1: FTIR Peak Assignments ---")
    print(f"The peak at 1645 cm^-1 corresponds to a {peak_assignments[1645]}.")
    print(f"The peak at 1652 cm^-1 corresponds to an {peak_assignments[1652]}.")
    print(f"The peaks at 1618 cm^-1 and 1680 cm^-1 are characteristic of a {peak_assignments[1618]}. The 1680 cm^-1 peak specifically suggests an anti-parallel arrangement.\n")

    # Step 2: Analyze the Concentration Titration Experiment
    # This experiment shows what happens upon gelation.
    print("--- Step 2: Analysis of Concentration Titration ---")
    print("The experiment shows a dual increase in the peaks at 1652 cm^-1 (Alpha Helix) and 1618 cm^-1 (Beta Sheet) as concentration increases.")
    print("This indicates that as the protein solution becomes more concentrated and forms a gel, the initially disordered proteins are folding to form both alpha-helical and beta-sheet structures.\n")

    # Step 3: Analyze the Heating Experiment
    # This experiment confirms the nature of the structures.
    print("--- Step 3: Analysis of Heating Experiment ---")
    print("Upon heating, the peaks at 1618 cm^-1 and 1680 cm^-1 (Beta Sheets) disappear.")
    print("Simultaneously, the peak at 1645 cm^-1 (Random Coil) grows stronger.")
    print("This shows that the ordered beta-sheet structures formed during gelation are thermally unstable and unfold back into a disordered state when heated. This process is called thermal denaturation.\n")
    
    # Step 4: Final Conclusion
    print("--- Step 4: Conclusion ---")
    print("The tardigrade proteins start in a disordered state. Upon hydration and an increase in concentration, they undergo a folding transition to form a hydrogel.")
    print("The FTIR data clearly shows the concurrent formation of both Alpha Helices (indicated by the 1652 cm^-1 peak) and Beta Sheets (indicated by the 1618 cm^-1 and 1680 cm^-1 peaks).")
    print("Therefore, the most accurate description of the gelation process is that disordered structures fold into both beta sheets and alpha helices.\n")

if __name__ == '__main__':
    analyze_ftir_data()
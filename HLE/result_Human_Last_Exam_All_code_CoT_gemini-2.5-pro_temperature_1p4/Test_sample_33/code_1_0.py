import sys

def analyze_ftir_data():
    """
    Analyzes FTIR data of tardigrade proteins to determine structural changes upon gelation.
    """
    # Define known FTIR peak assignments for protein secondary structures (Amide I band)
    peak_assignments = {
        "alpha-helix": "approx. 1650-1658 cm^-1",
        "beta-sheet (low frequency)": "approx. 1610-1640 cm^-1",
        "beta-sheet (high frequency, anti-parallel)": "approx. 1680-1695 cm^-1",
        "disordered/random coil": "approx. 1640-1650 cm^-1 (broad)"
    }

    # Experimental Observations from the problem
    observations = {
        "Concentration Titration (Gelation)": [
            "Dual increase in peak at 1652 cm^-1 as concentration increases.",
            "Dual increase in peak at 1618 cm^-1 as concentration increases."
        ],
        "Heating (De-gelation)": [
            "Peak at 1645 cm^-1 (broad) grows stronger.",
            "Peaks at 1618 cm^-1 and 1680 cm^-1 disappear."
        ]
    }

    print("Step 1: Assign protein structures to the observed FTIR peaks.")
    print(f"-> The peak at 1652 cm^-1 corresponds to an alpha-helix.")
    print(f"-> The peaks at 1618 cm^-1 and 1680 cm^-1 correspond to beta-sheets (specifically anti-parallel).")
    print(f"-> The peak at 1645 cm^-1 corresponds to a disordered (unfolded) structure.\n")

    print("Step 2: Analyze the concentration titration experiment, which shows structure formation.")
    print("The protein starts as disordered and forms a gel as concentration increases.")
    print(f"-> The observation that the 1652 cm^-1 peak (alpha-helix) increases means alpha-helices are forming.")
    print(f"-> The observation that the 1618 cm^-1 peak (beta-sheet) also increases means beta-sheets are forming.")
    print("Conclusion from this experiment: The gelation process involves the disordered protein folding into BOTH alpha-helices and beta-sheets.\n")
    
    print("Step 3: Analyze the heating experiment, which shows structure disruption.")
    print(f"-> The observation that the 1645 cm^-1 peak (disordered) increases confirms the protein is unfolding with heat.")
    print(f"-> The disappearance of the 1618 cm^-1 and 1680 cm^-1 peaks confirms the beta-sheet structures are being destroyed.")
    print("This confirms that the gel is made of ordered structures that unfold back to a disordered state upon heating.\n")

    print("Step 4: Final Conclusion.")
    print("The combined evidence shows that the initial disordered proteins fold to create a mixture of alpha-helical and beta-sheet structures to form the hydrogel.")
    
    # Redirect final answer to stdout
    sys.stdout.write("<<<I>>>")

if __name__ == '__main__':
    analyze_ftir_data()
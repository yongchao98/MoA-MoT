def analyze_enzyme_data():
    """
    Analyzes experimental data for enzyme Zma1 to determine the function
    of various molecules and selects the best-fitting conclusion.
    """
    experiments = {
        "Control": {"kcat": 500, "notes": "Zma1 + standard buffer"},
        "MgCl2": {"kcat": 700, "notes": "Control + 5 mM MgCl2"},
        "CaCl2": {"kcat": 500, "notes": "Control + 5 mM CaCl2"},
        "CuCl2": {"kcat": 400, "notes": "Control + 5 mM CuCl2"},
        "Al1": {"kcat": 1000, "notes": "Control + 5 mM Al1"},
        "Al2": {"kcat": 150, "notes": "Control + 5 mM Al2"},
        "Al1_Al2": {"kcat": 150, "notes": "Control + 5mM Al1 + 5mM Al2"},
        "XAG1": {"kcat": 10, "notes": "Control + 100 mM XAG1"},
        "XAG1_high_A": {"kcat": 450, "notes": "XAG1 condition + 500 mM of molecule A"},
        "Rga1": {"kcat": 10, "notes": "Control + 100 mM Rga1"},
        "Rga1_high_A": {"kcat": 10, "notes": "Rga1 condition + 500 mM of molecule A"}
    }

    control_kcat = experiments["Control"]["kcat"]

    print("--- Analysis of Experimental Data ---")

    # 1. Determine the function of Al1
    kcat_al1 = experiments["Al1"]["kcat"]
    print(f"\nFunction of Al1:")
    print(f"The kcat in the presence of Al1 is {kcat_al1}/s, compared to the control rate of {control_kcat}/s.")
    print("Conclusion: Al1 is an activator, likely an allosteric activator.")

    # 2. Determine the function of Rga1
    kcat_rga1 = experiments["Rga1"]["kcat"]
    kcat_rga1_high_A = experiments["Rga1_high_A"]["kcat"]
    kcat_xag1_high_A = experiments["XAG1_high_A"]["kcat"]
    print(f"\nFunction of Rga1:")
    print(f"Rga1 reduces kcat from {control_kcat}/s to {kcat_rga1}/s, so it is an inhibitor.")
    print(f"When substrate is increased, the kcat for Rga1 remains at {kcat_rga1_high_A}/s.")
    print(f"In contrast, for inhibitor XAG1, high substrate recovers the kcat to {kcat_xag1_high_A}/s.")
    print("Conclusion: The inability to reverse Rga1's inhibition with excess substrate strongly suggests it is a non-competitive or irreversible inhibitor.")

    # 3. Analyze the interaction of Al1 and Al2
    kcat_al2 = experiments["Al2"]["kcat"]
    kcat_al1_al2 = experiments["Al1_Al2"]["kcat"]
    print(f"\nInteraction of Al1 and Al2:")
    print(f"Al1 is an activator (kcat={kcat_al1}/s) and Al2 is an inhibitor (kcat={kcat_al2}/s).")
    print(f"When combined, their kcat is {kcat_al1_al2}/s, which is the same as Al2 alone.")
    print("Conclusion: This suggests they compete for the same binding site, and Al2's inhibitory effect is dominant.")

    print("\n--- Final Summary ---")
    print("Based on the data:")
    print("- Al1 is an allosteric activator.")
    print("- Rga1 is an irreversible inhibitor.")
    print("- Al1 and Al2 appear to compete for the same site.")

    print("\nEvaluating the choices, the one that best matches this analysis is C.")

if __name__ == '__main__':
    analyze_enzyme_data()
    print("\n<<<C>>>")

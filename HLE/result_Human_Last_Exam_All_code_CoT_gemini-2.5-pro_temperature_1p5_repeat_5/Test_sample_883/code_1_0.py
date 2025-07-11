def analyze_enzyme_data():
    """
    Analyzes experimental data for enzyme Zma1 to determine the function
    of various molecules and select the correct multiple-choice answer.
    """
    # Store the experimental results in a dictionary
    data = {
        "Control": 500,
        "MgCl2": 700,
        "CaCl2": 500,
        "CuCl2": 400,
        "Al1": 1000,
        "Al2": 150,
        "Al1_Al2": 150,
        "XAG1": 10,
        "XAG1_high_A": 450,
        "Rga1": 10,
        "Rga1_high_A": 10,
    }
    control_kcat = data["Control"]

    # --- Analysis for the primary question ---

    # 1. Analyze the function of Molecule Al1
    print("--- Analysis of Al1 ---")
    al1_kcat = data["Al1"]
    print(f"The control kcat is {control_kcat}/s.")
    print(f"In the presence of Al1, the kcat is {al1_kcat}/s.")
    print(f"Since {al1_kcat} > {control_kcat}, Al1 functions as an activator of Zma1.")
    print("An activator that binds to a site other than the active site is an allosteric regulator.\n")

    # 2. Analyze the function of Molecule Rga1
    print("--- Analysis of Rga1 ---")
    rga1_kcat = data["Rga1"]
    rga1_high_A_kcat = data["Rga1_high_A"]
    print(f"In the presence of Rga1, the kcat drops from {control_kcat}/s to {rga1_kcat}/s, indicating it is an inhibitor.")
    print(f"To determine the type of inhibition, we observe the effect of adding excess substrate (molecule A).")
    print(f"With Rga1 and high substrate, the kcat remains at {rga1_high_A_kcat}/s.")
    print(f"Since activity is not restored by excess substrate ({rga1_high_A_kcat}/s is not significantly higher than {rga1_kcat}/s), the inhibition is not competitive.")
    print("This is characteristic of non-competitive or irreversible inhibition. Among the choices, 'irreversible inhibitor' is the most fitting description.\n")

    # --- Evaluation of Answer Choice C ---
    print("--- Evaluation of Multiple Choice Options ---")
    print("Based on our analysis, let's evaluate choice C, which appears most promising.")
    print("Choice C: Al1 and Al2 function as allosteric modulators for Zma1. Al1 and Al2 bind to the same site on Zma1. Rga1 is an irreversible inhibitor.\n")

    # Part 1: Al1 and Al2 as allosteric modulators
    al2_kcat = data["Al2"]
    print("Verifying 'Al1 and Al2 function as allosteric modulators':")
    print(f"Al1 increases kcat from {control_kcat} to {al1_kcat} (activator).")
    print(f"Al2 decreases kcat from {control_kcat} to {al2_kcat} (inhibitor).")
    print("Conclusion: Correct. They both modulate the enzyme's activity.\n")

    # Part 2: Al1 and Al2 bind to the same site
    al1_al2_kcat = data["Al1_Al2"]
    print("Verifying 'Al1 and Al2 bind to the same site on Zma1':")
    print(f"With Al2 alone, kcat is {al2_kcat}/s.")
    print(f"When both Al1 (activator) and Al2 (inhibitor) are present, the kcat is {al1_al2_kcat}/s.")
    print(f"The resulting rate ({al1_al2_kcat}/s) is identical to the rate with the inhibitor Al2 alone ({al2_kcat}/s).")
    print("Conclusion: Correct. This suggests they compete for the same site, and the inhibitor's effect is dominant.\n")

    # Part 3: Rga1 is an irreversible inhibitor
    print("Verifying 'Rga1 is an irreversible inhibitor':")
    print("As determined in our initial analysis, the inhibition by Rga1 is not overcome by high substrate concentration.")
    print(f"(kcat of {rga1_kcat}/s with Rga1 vs. {rga1_high_A_kcat}/s with Rga1 + high substrate).")
    print("Conclusion: Correct. This matches the definition of irreversible or non-competitive inhibition.\n")

    print("--- Final Conclusion ---")
    print("All three statements in choice C are strongly supported by the experimental data.")


if __name__ == '__main__':
    analyze_enzyme_data()
<<<C>>>
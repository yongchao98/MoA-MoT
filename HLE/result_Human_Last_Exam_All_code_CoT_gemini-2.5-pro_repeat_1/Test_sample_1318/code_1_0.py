def analyze_phage_experiments():
    """
    Analyzes the provided experimental data and prints a step-by-step conclusion.
    """

    # --- Experiment 1 Data: Plaque-Forming Units (PFU) ---
    pfu_no_rp_wt = 100000
    pfu_no_rp_delta = 100000
    pfu_with_rp_wt = 80000
    pfu_with_rp_delta = 40000

    print("### Analysis of Experiment 1: Plaque-Forming Units (PFU) ###\n")

    # 1. Does RP system confer resistance?
    print("Step 1: Assessing the effect of the RP defense system.")
    print(f"- PFU for wild-type phage drops from {pfu_no_rp_wt} (no RP) to {pfu_with_rp_wt} (with RP).")
    print(f"- PFU for deltaXY phage drops from {pfu_no_rp_delta} (no RP) to {pfu_with_rp_delta} (with RP).")
    print("Conclusion: The RP system increases bacterial resistance to the phage.\n")

    # 2. What is the role of the XY operon?
    print("Step 2: Assessing the role of the phage's XY operon.")
    print(f"- In the presence of the RP system, the wild-type phage (with XY) forms {pfu_with_rp_wt} PFU.")
    print(f"- In the presence of the RP system, the deltaXY phage (without XY) forms {pfu_with_rp_delta} PFU.")
    improvement_factor = pfu_with_rp_wt / pfu_with_rp_delta
    print(f"Conclusion: The wild-type phage is {improvement_factor:.1f} times more effective. The XY operon helps the phage overcome the RP defense system.\n")

    # --- Experiment 2 Data: Mass Spectrometry ---
    # Using True for 'detected' and False for 'not detected'
    mass_spec = {
        "Sample 1 (RP + wt)": True,
        "Sample 2 (RP + deltaXY)": False,
        "Sample 3 (no RP + wt)": False,
        "Sample 4 (no RP + deltaXY)": False,
    }

    print("### Analysis of Experiment 2: Mass Spectrometry (500 Da molecule) ###\n")
    print("Step 3: Determining the origin of the 500 Da molecule.")
    print("Observation: The molecule is only detected after 60 minutes in Sample 1.")
    print("Sample 1 contains: Bacteria with RP system AND Phage with XY operon.")
    print("Conclusion: The production of the 500 Da molecule requires BOTH the host's RP system and the phage's XY operon.\n")

    # --- Synthesis and Final Conclusion ---
    print("### Final Synthesis ###\n")
    print("Combining results, we see that:")
    print("1. System RP is an anti-phage defense.")
    print("2. The phage's XY operon provides a counter-defense.")
    print("3. The counter-defense works by using a component of the RP system to produce a 500 Da molecule.")

    print("\nEvaluating Statement H: 'System RP increases the resistance of the bacteria against phageDE3 because the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP.'")
    print(f"- Part 1 is TRUE: Resistance increases (e.g., PFU drops from {pfu_no_rp_wt} to {pfu_with_rp_wt} for the wild-type phage).")
    print("- Part 2 is TRUE: The 500 Da product is only synthesized in the presence of both RP and XY.")
    print("This statement correctly connects the key findings from both experiments.")


analyze_phage_experiments()
<<<H>>>
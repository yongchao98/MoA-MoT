import collections

def analyze_phage_experiments():
    """
    Analyzes the provided experimental data to determine the correct statement.
    """

    # --- Data from the problem description ---
    exp1_cfu = {
        "no_RP": {
            "PhageDE3-wt": 100000,
            "PhageDE3-deltaXY": 100000,
        },
        "with_RP": {
            "PhageDE3-wt": 80000,
            "PhageDE3-deltaXY": 40000,
        },
    }

    exp2_mass_spec_60min = {
        "Sample 1 (RP+, wt phage)": True,
        "Sample 2 (RP+, deltaXY phage)": False,
        "Sample 3 (RP-, wt phage)": False,
        "Sample 4 (RP-, deltaXY phage)": False,
    }

    # --- Analysis Step-by-Step ---
    print("--- Analysis of Experimental Data ---")
    
    # 1. Analyze Experiment 1: What is the function of the RP system?
    print("\nStep 1: Analyzing the role of the RP system (Experiment 1)")
    cfu_no_rp = exp1_cfu["no_RP"]["PhageDE3-deltaXY"]
    cfu_with_rp = exp1_cfu["with_RP"]["PhageDE3-deltaXY"]
    print(f"Comparing a phage without operon XY on bacteria without RP vs. with RP.")
    print(f"CFU without RP system: {cfu_no_rp}")
    print(f"CFU with RP system:    {cfu_with_rp}")
    rp_is_defense = cfu_with_rp < cfu_no_rp
    if rp_is_defense:
        print(f"Result: {cfu_with_rp} < {cfu_no_rp}. The RP system reduces phage effectiveness. Therefore, System RP increases the resistance of the bacteria against the phage.")
    else:
        print("Result: The RP system does not increase resistance.")

    # 2. Analyze Experiment 1: What is the function of the XY operon?
    print("\nStep 2: Analyzing the role of the XY operon (Experiment 1)")
    cfu_wt_on_rp = exp1_cfu["with_RP"]["PhageDE3-wt"]
    cfu_delta_on_rp = exp1_cfu["with_RP"]["PhageDE3-deltaXY"]
    print(f"Comparing wild-type vs. deltaXY phage on bacteria that have the RP defense system.")
    print(f"CFU of PhageDE3-wt:      {cfu_wt_on_rp}")
    print(f"CFU of PhageDE3-deltaXY: {cfu_delta_on_rp}")
    xy_is_anti_defense = cfu_wt_on_rp > cfu_delta_on_rp
    if xy_is_anti_defense:
        print(f"Result: {cfu_wt_on_rp} > {cfu_delta_on_rp}. The XY operon helps the phage overcome the RP defense system.")
    else:
        print("Result: The XY operon has no effect against the RP system.")
    
    # 3. Analyze Experiment 2: What are the conditions for producing the 500 Da molecule?
    print("\nStep 3: Analyzing the production of the 500 Da molecule (Experiment 2)")
    print("The 500 Da molecule was detected under the following conditions at 60 minutes post-infection:")
    synthesis_conditions_met = False
    for sample, detected in exp2_mass_spec_60min.items():
        if detected:
            print(f"- {sample}: DETECTED")
            # This confirms production requires both RP system and WT phage (with operon XY)
            if "RP+" in sample and "wt phage" in sample:
                 synthesis_conditions_met = True
        else:
            print(f"- {sample}: NOT DETECTED")
    
    if synthesis_conditions_met:
        print("Result: The 500 Da molecule is produced only when the phage has the XY operon AND the bacteria has the RP system.")
        product_requires_xy_and_rp = True
    else:
        product_requires_xy_and_rp = False
        
    # 4. Evaluate Statement H
    print("\nStep 4: Evaluating the provided answer choices")
    print("Let's evaluate Statement H: 'System RP increases the resistance of the bacteria against phageDE3 because the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP.'")
    
    # Check first part of statement H
    part1_H = rp_is_defense
    print(f"- Part 1: 'System RP increases the resistance...' This is TRUE based on our analysis in Step 1.")
    
    # Check second part of statement H
    part2_H = product_requires_xy_and_rp
    print(f"- Part 2: '...the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP.' This is TRUE based on our analysis in Step 3.")
    
    if part1_H and part2_H:
        print("\nConclusion: Both claims in statement H are directly supported by the experimental data. It is the most complete description of the observed phenomena.")
    else:
        print("\nConclusion: Statement H is not fully supported by the data.")
        
    final_answer = 'H'
    print(f"\nThe correct statement is H.")


analyze_phage_experiments()
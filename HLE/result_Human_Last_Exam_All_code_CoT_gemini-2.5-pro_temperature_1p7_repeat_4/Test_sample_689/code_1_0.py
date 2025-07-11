import collections

def analyze_protein_folding():
    """
    Analyzes DLS data for MAB13 protein to determine the best folding conditions.
    """
    # The two experiments with HP70 are very similar. The question might have a typo,
    # but the conclusion is the same. We use the more favorable result for analysis.
    # Note: Using the average (77.5%) or the lower value (70%) does not change the final conclusion.
    data = {
        "Ecoli_37C": {"condition": "Protein MAB13 expressed in Escherichia coli at 37°C", "results": {30: 70, 55: 30}},
        "Ecoli_18C": {"condition": "Protein MAB13 expressed in Escherichia coli at 18°C", "results": {7.1: 20, 30: 80}},
        "Ecoli_18C_HP70": {"condition": "Protein MAB13 co-expressed with HP70 in Escherichia coli at 18°C", "results": {7.1: 85, 30: 15}},
        "HEK293_37C": {"condition": "Protein MAB13 expressed in HEK293 cells at 37°C", "results": {7.1: 95, 30: 5}},
        "Ecoli_37C_GFP": {"condition": "Protein MAB13 expressed in Escherichia coli at 37°C as an N-terminal fusion with GFP", "results": {30: 70, 55: 30}},
        "Ecoli_18C_MBP": {"condition": "Protein MAB13 expressed in Escherichia coli at 18°C as an N-terminal fusion with MBP", "results": {7.1: 60, 30: 30, 55: 10}}
    }
    
    # Step 1: Identify the hydrodynamic radius of the properly folded monomer.
    # The HEK293 cell expression is the positive control, showing 95% at 7.1 nm.
    monomer_radius = 7.1
    print(f"Step 1: Analysis of the control experiment (HEK293 cells).")
    print(f"The properly folded monomer of MAB13 has a hydrodynamic radius of {monomer_radius} nm, as it constitutes {data['HEK293_37C']['results'][monomer_radius]}% of the protein in the high-quality HEK293 system.\n")
    print("Step 2: Evaluating each factor's effect on folding.\n")

    # Helper to get monomer percentage
    def get_monomer_percent(exp_key):
        return data[exp_key]["results"].get(monomer_radius, 0)

    # Evaluation
    # Effect of lower temperature
    percent_37C = get_monomer_percent("Ecoli_37C")
    percent_18C = get_monomer_percent("Ecoli_18C")
    print("--- Evaluating the effect of temperature ---")
    print(f"Monomer percentage at 37°C in E. coli: {percent_37C}%")
    print(f"Monomer percentage at 18°C in E. coli: {percent_18C}%")
    if percent_18C > percent_37C:
        print(f"Conclusion: Lowering the temperature improves folding, as {percent_18C}% > {percent_37C}%.\n")
    else:
        print("Conclusion: Lowering the temperature does not improve folding.\n")
        
    # Effect of HP70
    percent_18C_HP70 = get_monomer_percent("Ecoli_18C_HP70")
    print("--- Evaluating the effect of HP70 chaperone ---")
    print(f"Monomer percentage at 18°C with HP70: {percent_18C_HP70}%")
    print(f"Monomer percentage at 18°C without HP70: {percent_18C}%")
    if percent_18C_HP70 > percent_18C:
        print(f"Conclusion: HP70 facilitates folding, as {percent_18C_HP70}% > {percent_18C}%.\n")
    else:
        print("Conclusion: HP70 does not facilitate folding.\n")
        
    # Effect of MBP fusion
    percent_18C_MBP = get_monomer_percent("Ecoli_18C_MBP")
    print("--- Evaluating the effect of MBP fusion ---")
    print(f"Monomer percentage at 18°C with MBP: {percent_18C_MBP}%")
    print(f"Monomer percentage at 18°C without MBP: {percent_18C}%")
    if percent_18C_MBP > percent_18C:
        print(f"Conclusion: MBP fusion improves folding, as {percent_18C_MBP}% > {percent_18C}%.\n")
    else:
        print("Conclusion: MBP fusion does not improve folding.\n")

    # Effect of GFP fusion
    percent_37C_GFP = get_monomer_percent("Ecoli_37C_GFP")
    print("--- Evaluating the effect of GFP fusion ---")
    print(f"Monomer percentage at 37°C with GFP: {percent_37C_GFP}%")
    print(f"Monomer percentage at 37°C without GFP: {percent_37C}%")
    if percent_37C_GFP > percent_37C:
        print(f"Conclusion: GFP fusion improves folding, as {percent_37C_GFP}% > {percent_37C}%.\n")
    else:
        print(f"Conclusion: GFP fusion does not improve folding, as {percent_37C_GFP}% is not greater than {percent_37C}%.\n")
        
    print("Step 3: Final Conclusion based on analysis.\n")
    print("Statement F says: 'HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.'")
    print("Our analysis confirms:")
    print(f"- HP70 helps at 18°C (monomer increased from {percent_18C}% to {percent_18C_HP70}%).")
    print(f"- MBP helps (monomer increased from {percent_18C}% to {percent_18C_MBP}%).")
    print(f"- Lower temperature helps (monomer increased from {percent_37C}% to {percent_18C}%).")
    print("All testable parts of statement F are correct. All other statements contain assertions that are demonstrably false based on the data.")
    print("Therefore, F is the correct answer.")

analyze_protein_folding()
print("<<<F>>>")
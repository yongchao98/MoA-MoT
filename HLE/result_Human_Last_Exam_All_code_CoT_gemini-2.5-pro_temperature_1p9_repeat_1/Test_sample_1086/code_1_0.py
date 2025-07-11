import pandas as pd

# Store the experimental data in a structured DataFrame
data = {
    'Mouse Line': ['Wild-type', 'delta-ber1', 'delta-ber2', 'delta-ber1, delta-ber2'],
    'Center Time (%)': [15, 15, 8, 8],
    'Distance (cm)': [900, 900, 1250, 1250],
    'Immobility (s)': [180, 180, 230, 230],
    'Sucrose Pref (%)': [75, 62, 62, 62],
    'Ki67 Cells': [3500, 3500, 3500, 2850]
}
df = pd.DataFrame(data).set_index('Mouse Line')

# Store the SSRI experiment data
ssri_data = {
    'Mouse Line': ['Wild-type', 'delta-ber1', 'delta-ber2', 'delta-ber1, delta-ber2'],
    'Center Time (%)': [15, 15, 15, 15],
    'Distance (cm)': [900, 900, 900, 900]
}
ssri_df = pd.DataFrame(ssri_data).set_index('Mouse Line')

# Retrieve baseline values from Wild-type
wt_baseline = df.loc['Wild-type']
ssri_baseline_center_time = ssri_df.loc['Wild-type', 'Center Time (%)']
ssri_baseline_distance = ssri_df.loc['Wild-type', 'Distance (cm)']


print("Analyzing the statements in Answer A:")

print("\n1. 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI.'")
# Check the anxiety/hyperactivity phenotype in ber2 mutants
mutant_line = 'delta-ber2'
original_center_time = df.loc[mutant_line, 'Center Time (%)']
ssri_center_time = ssri_df.loc[mutant_line, 'Center Time (%)']
print(f"   - For {mutant_line}, SSRI treatment reversed the anxiety phenotype (Center Time): {original_center_time}% -> {ssri_center_time}%.")

original_distance = df.loc[mutant_line, 'Distance (cm)']
ssri_distance = ssri_df.loc[mutant_line, 'Distance (cm)']
print(f"   - For {mutant_line}, SSRI treatment reversed the hyperactivity phenotype (Distance): {original_distance} cm -> {ssri_distance} cm.")
print("   - This statement is supported by the data.")


print("\n2. 'Mice with defects in ber2 may not have a decrease in cell proliferation.'")
ber2_ki67 = df.loc['delta-ber2', 'Ki67 Cells']
wt_ki67 = wt_baseline['Ki67 Cells']
print(f"   - Ki67 cell count in delta-ber2 mice is {ber2_ki67}, which is the same as the Wild-type count of {wt_ki67}.")
print("   - This statement is true, as a ber2 defect alone does not cause a decrease.")


print("\n3. 'Gene ber1 and ber2 regulate cell proliferation.'")
ber1_ki67 = df.loc['delta-ber1', 'Ki67 Cells']
double_ko_ki67 = df.loc['delta-ber1, delta-ber2', 'Ki67 Cells']
print(f"   - Ki67 in delta-ber1 is {ber1_ki67} (no change from Wild-type).")
print(f"   - Ki67 in delta-ber2 is {ber2_ki67} (no change from Wild-type).")
print(f"   - However, Ki67 in the double knockout ('delta-ber1, delta-ber2') is {double_ko_ki67}, a decrease from the Wild-type count of {wt_ki67}.")
print("   - This indicates a redundant function, and that both genes are involved in regulating cell proliferation. The statement is true.")


print("\nBased on the analysis, all statements in answer A are correct and supported by the experimental data.")
print('<<<A>>>')
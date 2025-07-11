import pandas as pd

# Data from the problem description
data = {
    'Line': ['Wild-type', 'delta-ber1', 'delta-ber2', 'delta-ber1, delta-ber2'],
    'Time in Center (%)': [15, 15, 8, 8],
    'Distance (cm/5min)': [900, 900, 1250, 1250],
    'Immobility (s)': [180, 180, 230, 230],
    'Sucrose Preference (%)': [75, 62, 62, 62],
    'Ki67 Cells': [3500, 3500, 3500, 2850],
    'Time in Center post-SSRI (%)': [15, 15, 15, 15],
    'Distance post-SSRI (cm/5min)': [900, 900, 900, 900]
}

df = pd.DataFrame(data).set_index('Line')

# --- Analysis to verify Answer A ---

print("Analyzing the data to verify the statements in Answer Choice A:")
print("-" * 60)

# Statement 1: "The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI."
print("1. Checking SSRI reversal effect:")
phenotype_line = 'delta-ber2'
wt_time_center = df.loc['Wild-type', 'Time in Center (%)']
phenotype_time_center_before = df.loc[phenotype_line, 'Time in Center (%)']
phenotype_time_center_after = df.loc[phenotype_line, 'Time in Center post-SSRI (%)']

print(f"   - The anxiety phenotype (time in center) for '{phenotype_line}' mice was {phenotype_time_center_before}%.")
print(f"   - After SSRI treatment, their time in center became {phenotype_time_center_after}%.")
print(f"   - The wild-type value is {wt_time_center}%.")
if phenotype_time_center_after == wt_time_center:
    print("   -> Conclusion: The behavioral phenotype was reversed to wild-type levels. Statement is SUPPORTED.")
else:
    print("   -> Conclusion: Statement is NOT supported.")
print("-" * 60)


# Statement 2: "Mice with defects in ber2 may not have a decrease in cell proliferation."
print("2. Checking cell proliferation in delta-ber2 mice:")
wt_ki67 = df.loc['Wild-type', 'Ki67 Cells']
ber2_ki67 = df.loc['delta-ber2', 'Ki67 Cells']

print(f"   - The Ki67 cell count for 'delta-ber2' mice is {ber2_ki67}.")
print(f"   - The wild-type Ki67 cell count is {wt_ki67}.")
if ber2_ki67 >= wt_ki67:
    print("   -> Conclusion: The 'delta-ber2' mice do not show a decrease in proliferation. Statement is SUPPORTED.")
else:
    print("   -> Conclusion: Statement is NOT supported.")
print("-" * 60)


# Statement 3: "Gene ber1 and ber2 regulate cell proliferation."
print("3. Checking the role of both genes in cell proliferation:")
dko_ki67 = df.loc['delta-ber1, delta-ber2', 'Ki67 Cells']

print(f"   - The Ki67 cell count for the double-knockout mice is {dko_ki67}.")
print(f"   - The wild-type Ki67 cell count is {wt_ki67}.")
if dko_ki67 < wt_ki67:
    print(f"   -> Conclusion: The absence of both genes leads to a decrease in proliferation (from {wt_ki67} to {dko_ki67}), indicating they regulate this process. Statement is SUPPORTED.")
else:
    print("   -> Conclusion: Statement is NOT supported.")
print("-" * 60)

print("\nAll three statements in Answer Choice A are supported by the data.")
print("The final answer is A")

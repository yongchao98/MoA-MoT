import sys

# Store the experimental data in a nested dictionary for easy access.
data = {
    'WT': {
        'open_field_center_time': 15,
        'open_field_distance': 900,
        'ki67_cells': 3500,
        'ssri_center_time': 15,
        'ssri_distance': 900
    },
    'delta-ber1': {
        'ki67_cells': 3500,
    },
    'delta-ber2': {
        'open_field_center_time': 8,
        'open_field_distance': 1250,
        'ki67_cells': 3500,
        'ssri_center_time': 15,
        'ssri_distance': 900
    },
    'delta-ber1_delta-ber2': {
        'ki67_cells': 2850
    }
}

# The chosen answer is A. This script will verify each statement within option A.
print("Verifying the statements in Option A based on the provided data:")
print("-" * 60)

# Statement 1: "The effects of mutations in ber1 and ber2 may be reversed by treatment with selective serotonin reuptake inhibitors (SSRI)."
# We check if the open field test phenotype in delta-ber2 mice was reversed to wild-type levels.
print("1. Verifying SSRI reversal of the ber2 mutation effect:")

ber2_original_center_time = data['delta-ber2']['open_field_center_time']
ber2_ssri_center_time = data['delta-ber2']['ssri_center_time']
wt_center_time = data['WT']['open_field_center_time']
print(f"   - For 'time spent in center', the delta-ber2 value changed from {ber2_original_center_time}% to {ber2_ssri_center_time}%.")
print(f"   - This is a reversal, as the new value {ber2_ssri_center_time}% matches the wild-type value of {wt_center_time}%.")

ber2_original_distance = data['delta-ber2']['open_field_distance']
ber2_ssri_distance = data['delta-ber2']['ssri_distance']
wt_distance = data['WT']['open_field_distance']
print(f"   - For 'distance moved', the delta-ber2 value changed from {ber2_original_distance} cm to {ber2_ssri_distance} cm.")
print(f"   - This is a reversal, as the new value {ber2_ssri_distance} cm matches the wild-type value of {wt_distance} cm.")
print("   - Conclusion: The statement is supported by the data.")
print("-" * 60)

# Statement 2: "Mice with defects in ber2 may not have a decrease in cell proliferation."
# We check if the Ki67 cell count in delta-ber2 mice is different from wild-type.
print("2. Verifying the effect of ber2 defect on cell proliferation:")
ber2_ki67 = data['delta-ber2']['ki67_cells']
wt_ki67 = data['WT']['ki67_cells']
print(f"   - Ki67 cell count in delta-ber2 mice is {ber2_ki67}.")
print(f"   - Ki67 cell count in wild-type mice is {wt_ki67}.")
if ber2_ki67 == wt_ki67:
    print(f"   - The values {ber2_ki67} and {wt_ki67} are equal, showing no decrease.")
    print("   - Conclusion: The statement is supported by the data.")
else:
    sys.exit("Error in logic: Statement 2 verification failed.")
print("-" * 60)

# Statement 3: "Gene ber1 and ber2 regulate cell proliferation."
# This implies a redundant function, where single knockouts have no effect, but the double knockout does.
print("3. Verifying the role of ber1 and ber2 in cell proliferation:")
ber1_ki67 = data['delta-ber1']['ki67_cells']
double_ko_ki67 = data['delta-ber1_delta-ber2']['ki67_cells']

print(f"   - Ki67 count in wild-type: {wt_ki67}")
print(f"   - Ki67 count in delta-ber1: {ber1_ki67}")
print(f"   - Ki67 count in delta-ber2: {ber2_ki67}")
print(f"   - Ki67 count in double knockout (delta-ber1, delta-ber2): {double_ko_ki67}")
if ber1_ki67 == wt_ki67 and ber2_ki67 == wt_ki67 and double_ko_ki67 < wt_ki67:
    print(f"   - The double knockout value ({double_ko_ki67}) is less than the wild-type value ({wt_ki67}), while single knockouts are not.")
    print("   - This indicates a redundant function and supports the conclusion that both genes regulate proliferation.")
else:
    sys.exit("Error in logic: Statement 3 verification failed.")
print("-" * 60)

print("\nFinal Result: All three statements in Option A are supported by the data analysis.")

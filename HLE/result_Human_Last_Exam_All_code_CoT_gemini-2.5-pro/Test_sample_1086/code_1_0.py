import collections

# Storing the experimental data in a structured way
data = {
    "Wild-type": {
        "open_field_center_time_percent": 15,
        "forced_swim_immobility_sec": 180,
        "sucrose_preference_percent": 75,
        "ki67_cells": 3500,
        "open_field_center_time_post_ssri_percent": 15
    },
    "delta-ber1": {
        "open_field_center_time_percent": 15,
        "forced_swim_immobility_sec": 180,
        "sucrose_preference_percent": 62,
        "ki67_cells": 3500,
        "open_field_center_time_post_ssri_percent": 15
    },
    "delta-ber2": {
        "open_field_center_time_percent": 8,
        "forced_swim_immobility_sec": 230,
        "sucrose_preference_percent": 62,
        "ki67_cells": 3500,
        "open_field_center_time_post_ssri_percent": 15
    },
    "delta-ber1, delta-ber2": {
        "open_field_center_time_percent": 8,
        "forced_swim_immobility_sec": 230,
        "sucrose_preference_percent": 62,
        "ki67_cells": 2850,
        "open_field_center_time_post_ssri_percent": 15
    }
}

# We will evaluate the statements in Answer Choice A.

print("Evaluating Answer Choice A based on the experimental data:")
print("-" * 50)

# Statement 1: "The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI."
print("1. Checking SSRI treatment effects:")
ber2_pre_ssri = data["delta-ber2"]["open_field_center_time_percent"]
ber2_post_ssri = data["delta-ber2"]["open_field_center_time_post_ssri_percent"]
double_ko_pre_ssri = data["delta-ber1, delta-ber2"]["open_field_center_time_percent"]
double_ko_post_ssri = data["delta-ber1, delta-ber2"]["open_field_center_time_post_ssri_percent"]
wt_value = data["Wild-type"]["open_field_center_time_percent"]

print(f"   - The anxiety phenotype (time in center) for delta-ber2 mice was {ber2_pre_ssri}% before SSRI treatment.")
print(f"   - After SSRI treatment, it became {ber2_post_ssri}%, matching the Wild-type value of {wt_value}%.")
print(f"   - Similarly, for the double knockout, the value changed from {double_ko_pre_ssri}% to {double_ko_post_ssri}%.")
print("   - Conclusion: This supports the statement that effects may be reversed.\n")


# Statement 2: "Mice with defects in ber2 may not have a decrease in cell proliferation."
print("2. Checking the effect of ber2 defect on cell proliferation:")
wt_ki67 = data["Wild-type"]["ki67_cells"]
ber2_ki67 = data["delta-ber2"]["ki67_cells"]

print(f"   - Wild-type mice had {wt_ki67} Ki67-positive cells.")
print(f"   - delta-ber2 mice had {ber2_ki67} Ki67-positive cells.")
print("   - Conclusion: The number of cells did not decrease. This supports the statement.\n")


# Statement 3: "Gene ber1 and ber2 regulate cell proliferation."
print("3. Checking the combined role of ber1 and ber2 on cell proliferation:")
ber1_ki67 = data["delta-ber1"]["ki67_cells"]
double_ko_ki67 = data["delta-ber1, delta-ber2"]["ki67_cells"]

print(f"   - Single knockouts of ber1 ({ber1_ki67} cells) and ber2 ({ber2_ki67} cells) did not change cell proliferation compared to wild-type ({wt_ki67} cells).")
print(f"   - However, the double knockout of ber1 and ber2 resulted in a decreased cell count of {double_ko_ki67}.")
print("   - Conclusion: This indicates a redundant or combined function, supporting that together they regulate cell proliferation.\n")


print("-" * 50)
print("All three statements in Answer Choice A are supported by the provided data.")

print("<<<A>>>")
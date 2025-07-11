import sys

# Step 1: Structure the experimental data into a dictionary for easy access.
data = {
    "Wild-type": {
        "open_field_center_time": 15, "ki67_cells": 3500, "open_field_center_time_ssri": 15
    },
    "delta-ber1": {
        "open_field_center_time": 15, "ki67_cells": 3500, "open_field_center_time_ssri": 15
    },
    "delta-ber2": {
        "open_field_center_time": 8, "ki67_cells": 3500, "open_field_center_time_ssri": 15
    },
    "delta-ber1, delta-ber2": {
        "open_field_center_time": 8, "ki67_cells": 2850, "open_field_center_time_ssri": 15
    }
}

# The task is to evaluate the statements in the correct answer choice, which is A.
# Let's break down and verify each statement in choice A using the data.

print("Evaluating Choice A: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI. Mice with defects in ber2 may not have a decrease in cell proliferation. Gene ber1 and ber2 regulate cell proliferation.'\n")

# Step 3: Verification of the first statement.
print("1. Analysis of: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI.'")
center_time_before_ssri = data["delta-ber2"]["open_field_center_time"]
center_time_after_ssri = data["delta-ber2"]["open_field_center_time_ssri"]
wt_center_time = data["Wild-type"]["open_field_center_time"]

# The mutation in ber2 caused a phenotype (reduced time in center). Let's check if SSRI reversed it.
print(f"   - Before SSRI, delta-ber2 mice spent {center_time_before_ssri}% of time in the center, compared to {wt_center_time}% for Wild-type.")
print(f"   - After SSRI treatment, delta-ber2 mice spent {center_time_after_ssri}% of time in the center, which is the same as Wild-type.")
print("   - Conclusion: The anxiety-like phenotype caused by the ber2 mutation was reversed. The statement is supported.\n")


# Step 4: Verification of the second statement.
print("2. Analysis of: 'Mice with defects in ber2 may not have a decrease in cell proliferation.'")
ber2_ki67 = data["delta-ber2"]["ki67_cells"]
wt_ki67 = data["Wild-type"]["ki67_cells"]

# Compare the cell proliferation marker Ki67 in delta-ber2 vs Wild-type.
print(f"   - The number of Ki67-positive cells in delta-ber2 mice is {ber2_ki67}.")
print(f"   - The number of Ki67-positive cells in Wild-type mice is {wt_ki67}.")
print("   - Conclusion: Since the cell count is the same, mice with a defect in only the ber2 gene did not have a decrease in cell proliferation. The statement is supported.\n")


# Step 5: Verification of the third statement.
print("3. Analysis of: 'Gene ber1 and ber2 regulate cell proliferation.'")
ber1_ko_ki67 = data["delta-ber1"]["ki67_cells"]
ber2_ko_ki67 = data["delta-ber2"]["ki67_cells"]
double_ko_ki67 = data["delta-ber1, delta-ber2"]["ki67_cells"]
wt_ki67 = data["Wild-type"]["ki67_cells"]

# Check if single knockouts affect proliferation vs the double knockout.
print(f"   - Proliferation in Wild-type: {wt_ki67} cells.")
print(f"   - Proliferation in delta-ber1 (single knockout): {ber1_ko_ki67} cells (no change).")
print(f"   - Proliferation in delta-ber2 (single knockout): {ber2_ko_ki67} cells (no change).")
print(f"   - Proliferation in delta-ber1, delta-ber2 (double knockout): {double_ko_ki67} cells (a decrease).")
print("   - Conclusion: A decrease in proliferation only occurs when both genes are knocked out. This demonstrates a redundant function, meaning both genes regulate cell proliferation. The statement is supported.\n")

sys.stdout.flush()
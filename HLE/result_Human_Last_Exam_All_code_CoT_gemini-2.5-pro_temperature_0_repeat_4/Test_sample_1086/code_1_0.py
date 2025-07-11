import sys

# Redirect print to a string to capture all output for the final response
class OutputCollector:
    def __init__(self):
        self.output = ""
    def write(self, text):
        self.output += text
    def flush(self):
        pass

original_stdout = sys.stdout
collector = OutputCollector()
sys.stdout = collector

# --- Data from the experiment ---
data = {
    "Wild-type": {
        "open_field_center": 15, "ki67_cells": 3500, "ssri_open_field_center": 15
    },
    "delta-ber1": {
        "ki67_cells": 3500
    },
    "delta-ber2": {
        "open_field_center": 8, "ki67_cells": 3500, "ssri_open_field_center": 15
    },
    "delta-ber1, delta-ber2": {
        "open_field_center": 8, "ki67_cells": 2850, "ssri_open_field_center": 15
    }
}

# --- Analysis of Statement A ---
print("Analyzing Statement A: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI. Mice with defects in ber2 may not have a decrease in cell proliferation. Gene ber1 and ber2 regulate cell proliferation.'")

# 1. Check Clause 1: SSRI Reversal
print("\n--- Step 1: Analyzing SSRI Treatment Effect ---")
wt_center_time = data["Wild-type"]["ssri_open_field_center"]
d_ber2_before = data["delta-ber2"]["open_field_center"]
d_ber2_after = data["delta-ber2"]["ssri_open_field_center"]
d_double_before = data["delta-ber1, delta-ber2"]["open_field_center"]
d_double_after = data["delta-ber1, delta-ber2"]["ssri_open_field_center"]

print(f"The anxiety phenotype in delta-ber2 mice (time in center) was {d_ber2_before}%. After SSRI treatment, it became {d_ber2_after}%.")
print(f"The anxiety phenotype in delta-ber1, delta-ber2 mice was {d_double_before}%. After SSRI treatment, it became {d_double_after}%.")
print(f"The wild-type level is {wt_center_time}%. Since the values in mutant mice returned to wild-type levels, the phenotype was reversed.")
print("Conclusion: The first clause is supported.")

# 2. Check Clause 2: Cell Proliferation in delta-ber2
print("\n--- Step 2: Analyzing Cell Proliferation in delta-ber2 Mice ---")
wt_ki67 = data["Wild-type"]["ki67_cells"]
d_ber2_ki67 = data["delta-ber2"]["ki67_cells"]

print(f"Wild-type mice have {wt_ki67} Ki67+ cells.")
print(f"delta-ber2 mice have {d_ber2_ki67} Ki67+ cells.")
print(f"Since {d_ber2_ki67} is equal to {wt_ki67}, a defect in ber2 alone does not cause a decrease in cell proliferation.")
print("Conclusion: The second clause is supported.")

# 3. Check Clause 3: Redundant Function in Proliferation
print("\n--- Step 3: Analyzing the Role of ber1 and ber2 in Proliferation ---")
d_ber1_ki67 = data["delta-ber1"]["ki67_cells"]
d_double_ki67 = data["delta-ber1, delta-ber2"]["ki67_cells"]

print(f"Ki67+ cells count:")
print(f"  - Wild-type: {wt_ki67}")
print(f"  - delta-ber1 (single knockout): {d_ber1_ki67}")
print(f"  - delta-ber2 (single knockout): {d_ber2_ki67}")
print(f"  - delta-ber1, delta-ber2 (double knockout): {d_double_ki67}")
print(f"Proliferation is normal in single knockouts but is reduced from {wt_ki67} to {d_double_ki67} when both genes are knocked out.")
print("This indicates a combined or redundant function of ber1 and ber2 in regulating cell proliferation.")
print("Conclusion: The third clause is supported.")

# 4. Final Conclusion
print("\n--- Final Conclusion ---")
print("All three clauses of statement A are strongly supported by the experimental data.")

# Restore original stdout and print the collected output
sys.stdout = original_stdout
print(collector.output)
print("<<<A>>>")
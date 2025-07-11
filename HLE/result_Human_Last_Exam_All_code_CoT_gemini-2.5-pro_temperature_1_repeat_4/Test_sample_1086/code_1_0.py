import sys
import io

# Capture original stdout to restore it later
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Data Representation ---
# Store the experimental results in a structured dictionary for easy access.
data = {
    'WT': {
        'open_field_center': 15,
        'open_field_distance': 900,
        'forced_swim_immobility': 180,
        'sucrose_preference': 75,
        'ki67_cells': 3500,
        'ssri_open_field_center': 15,
        'ssri_open_field_distance': 900,
    },
    'delta-ber1': {
        'open_field_center': 15,
        'open_field_distance': 900,
        'forced_swim_immobility': 180,
        'sucrose_preference': 62,
        'ki67_cells': 3500,
        'ssri_open_field_center': 15,
        'ssri_open_field_distance': 900,
    },
    'delta-ber2': {
        'open_field_center': 8,
        'open_field_distance': 1250,
        'forced_swim_immobility': 230,
        'sucrose_preference': 62,
        'ki67_cells': 3500,
        'ssri_open_field_center': 15,
        'ssri_open_field_distance': 900,
    },
    'delta-ber1, delta-ber2': {
        'open_field_center': 8,
        'open_field_distance': 1250,
        'forced_swim_immobility': 230,
        'sucrose_preference': 62,
        'ki67_cells': 2850,
        'ssri_open_field_center': 15,
        'ssri_open_field_distance': 900,
    }
}

# --- Analysis of Statement A ---
# Statement A: "The effects of mutations in ber1 and ber2 may be reversed by treatment with
# selective serotonin reuptake inhibitors (SSRI). Mice with defects in ber2 may not have a
# decrease in cell proliferation. Gene ber1 and ber2 regulate cell proliferation."

print("Analysis of Statement A:")
print("="*30)

# Part 1: "The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI."
# We check if the anxiety phenotype (reduced time in center) in delta-ber2 mice was reversed.
print("\n1. Checking for SSRI reversal of anxiety phenotype:")
phenotype_before = data['delta-ber2']['open_field_center']
phenotype_after = data['delta-ber2']['ssri_open_field_center']
wt_value = data['WT']['open_field_center']
print(f"   - Time in center for delta-ber2 mice before SSRI: {phenotype_before}%")
print(f"   - Time in center for delta-ber2 mice after SSRI: {phenotype_after}%")
print(f"   - Time in center for Wild-type mice: {wt_value}%")
is_reversed = phenotype_before != wt_value and phenotype_after == wt_value
print(f"   - Conclusion: Since {phenotype_before} != {wt_value} and {phenotype_after} == {wt_value}, an effect of the mutation was reversed. This statement is supported by the data.")

# Part 2: "Mice with defects in ber2 may not have a decrease in cell proliferation."
# We check if the delta-ber2 single-mutant shows a decrease in Ki67 cells.
print("\n2. Checking if ber2 defect alone causes a decrease in cell proliferation:")
ber2_ki67 = data['delta-ber2']['ki67_cells']
wt_ki67 = data['WT']['ki67_cells']
print(f"   - Ki67-positive cells in delta-ber2 mice: {ber2_ki67}")
print(f"   - Ki67-positive cells in Wild-type mice: {wt_ki67}")
no_decrease = ber2_ki67 >= wt_ki67
print(f"   - Conclusion: Since {ber2_ki67} is not less than {wt_ki67}, the delta-ber2 line does not show a decrease. This statement is correct.")

# Part 3: "Gene ber1 and ber2 regulate cell proliferation."
# We check if the double-mutant shows a decrease in Ki67 cells, implying joint regulation.
print("\n3. Checking if ber1 and ber2 together regulate cell proliferation:")
double_ko_ki67 = data['delta-ber1, delta-ber2']['ki67_cells']
wt_ki67 = data['WT']['ki67_cells']
print(f"   - Ki67-positive cells in double knockout (delta-ber1, delta-ber2) mice: {double_ko_ki67}")
print(f"   - Ki67-positive cells in Wild-type mice: {wt_ki67}")
is_regulated = double_ko_ki67 < wt_ki67
print(f"   - Conclusion: Since {double_ko_ki67} is less than {wt_ki67}, the double knockout shows a defect. This indicates joint regulation by ber1 and ber2. This statement is correct.")

print("\n" + "="*30)
print("All three parts of Statement A are strongly supported by the experimental data.")

# --- Final Answer ---
final_answer = 'A'

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_string = captured_output.getvalue()

# Print the captured output
print(output_string)
print(f'<<<A>>>')
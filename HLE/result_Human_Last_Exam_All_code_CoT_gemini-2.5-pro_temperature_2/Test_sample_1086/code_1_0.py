# Experimental data stored in a dictionary
data = {
    'Wild-type': {
        'open_field_center_time_pre': 15, 'open_field_distance_pre': 900,
        'immobility_time': 180, 'sucrose_preference': 75, 'ki67_cells': 3500,
        'open_field_center_time_post_ssri': 15, 'open_field_distance_post_ssri': 900,
    },
    'delta-ber1': {
        'open_field_center_time_pre': 15, 'open_field_distance_pre': 900,
        'immobility_time': 180, 'sucrose_preference': 62, 'ki67_cells': 3500,
        'open_field_center_time_post_ssri': 15, 'open_field_distance_post_ssri': 900,
    },
    'delta-ber2': {
        'open_field_center_time_pre': 8, 'open_field_distance_pre': 1250,
        'immobility_time': 230, 'sucrose_preference': 62, 'ki67_cells': 3500,
        'open_field_center_time_post_ssri': 15, 'open_field_distance_post_ssri': 900,
    },
    'delta-ber1, delta-ber2': {
        'open_field_center_time_pre': 8, 'open_field_distance_pre': 1250,
        'immobility_time': 230, 'sucrose_preference': 62, 'ki67_cells': 2850,
        'open_field_center_time_post_ssri': 15, 'open_field_distance_post_ssri': 900,
    }
}

# Aliases for easier access
wt = data['Wild-type']
d_ber2 = data['delta-ber2']
d_ber1_ber2 = data['delta-ber1, delta-ber2']

print("Analyzing the statements in answer choice A:\n")

# --- Statement 1 Analysis ---
print("1. 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI.'")
print("   - In delta-ber2 mice, the time in the center was {}% before SSRI and {}% after SSRI (WT is {}%).".format(d_ber2['open_field_center_time_pre'], d_ber2['open_field_center_time_post_ssri'], wt['open_field_center_time_pre']))
print("   - The distance moved was {} cm before and {} cm after (WT is {} cm).".format(d_ber2['open_field_distance_pre'], d_ber2['open_field_distance_post_ssri'], wt['open_field_distance_pre']))
print("   - In delta-ber1, delta-ber2 mice, the time in the center was {}% before SSRI and {}% after SSRI.".format(d_ber1_ber2['open_field_center_time_pre'], d_ber1_ber2['open_field_center_time_post_ssri']))
print("   - The distance moved was {} cm before and {} cm after.".format(d_ber1_ber2['open_field_distance_pre'], d_ber1_ber2['open_field_distance_post_ssri']))
if (d_ber2['open_field_center_time_post_ssri'] == wt['open_field_center_time_pre'] and
    d_ber1_ber2['open_field_center_time_post_ssri'] == wt['open_field_center_time_pre']):
    print("   -> Conclusion: The anxiety and hyperactivity phenotype in mice with a ber2 mutation was reversed. This statement is supported by the data.\n")
else:
    print("   -> Conclusion: This statement is not supported.\n")


# --- Statement 2 Analysis ---
print("2. 'Mice with defects in ber2 may not have a decrease in cell proliferation.'")
print("   - Wild-type mice have {} Ki67 cells.".format(wt['ki67_cells']))
print("   - delta-ber2 mice have {} Ki67 cells.".format(d_ber2['ki67_cells']))
if d_ber2['ki67_cells'] == wt['ki67_cells']:
    print("   -> Conclusion: A defect in only the ber2 gene did not cause a decrease in cell proliferation. This statement is correct.\n")
else:
    print("   -> Conclusion: This statement is incorrect.\n")

# --- Statement 3 Analysis ---
print("3. 'Gene ber1 and ber2 regulate cell proliferation.'")
print("   - Wild-type mice have {} Ki67 cells.".format(wt['ki67_cells']))
print("   - delta-ber1, delta-ber2 (double knockout) mice have {} Ki67 cells.".format(d_ber1_ber2['ki67_cells']))
if d_ber1_ber2['ki67_cells'] < wt['ki67_cells']:
     print("   -> Conclusion: Since knocking out both genes causes a decrease in Ki67 cells, it demonstrates they are involved in regulating cell proliferation. This statement is correct.")
else:
     print("   -> Conclusion: This statement is incorrect.")

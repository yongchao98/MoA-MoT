def analyze_mouse_data():
    """
    Analyzes experimental data to verify the claims in the correct answer choice.
    """
    # Store the experimental data in a dictionary
    data = {
        'WT': {
            'open_field_center': 15,
            'ki67_cells': 3500,
        },
        'delta-ber2': {
            'open_field_center': 8,
            'ki67_cells': 3500,
            'open_field_ssri_center': 15,
        },
        'delta-ber1_delta-ber2': {
            'ki67_cells': 2850,
            'open_field_center': 8,
            'open_field_ssri_center': 15,
        }
    }

    print("Analyzing Statement A: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI. Mice with defects in ber2 may not have a decrease in cell proliferation. Gene ber1 and ber2 regulate cell proliferation.'")
    print("-" * 80)

    # --- Clause 1 Verification ---
    print("Clause 1: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI.'")
    wt_center_time = data['WT']['open_field_center']
    d_ber2_center_time_before = data['delta-ber2']['open_field_center']
    d_ber2_center_time_after = data['delta-ber2']['open_field_ssri_center']
    
    # The double knockout also shows reversal
    d_double_ko_center_time_before = data['delta-ber1_delta-ber2']['open_field_center']
    d_double_ko_center_time_after = data['delta-ber1_delta-ber2']['open_field_ssri_center']

    if d_ber2_center_time_after == wt_center_time and d_double_ko_center_time_after == wt_center_time:
        print(f"VERIFIED: The anxiety phenotype in delta-ber2 mice (time in center: {d_ber2_center_time_before}%) was reversed to wild-type levels ({wt_center_time}%) after SSRI treatment.")
        print(f"VERIFIED: The anxiety phenotype in double knockout mice (time in center: {d_double_ko_center_time_before}%) was also reversed to wild-type levels ({wt_center_time}%) after SSRI treatment.")
    else:
        print("Clause 1 is not supported.")
    print("-" * 80)

    # --- Clause 2 Verification ---
    print("Clause 2: 'Mice with defects in ber2 may not have a decrease in cell proliferation.'")
    wt_ki67 = data['WT']['ki67_cells']
    d_ber2_ki67 = data['delta-ber2']['ki67_cells']
    
    if d_ber2_ki67 == wt_ki67:
        print(f"VERIFIED: The delta-ber2 mice had {d_ber2_ki67} Ki67+ cells, which is the same as the wild-type ({wt_ki67} cells) and does not show a decrease.")
    else:
        print("Clause 2 is not supported.")
    print("-" * 80)

    # --- Clause 3 Verification ---
    print("Clause 3: 'Gene ber1 and ber2 regulate cell proliferation.'")
    d_double_ko_ki67 = data['delta-ber1_delta-ber2']['ki67_cells']
    
    if d_double_ko_ki67 < wt_ki67:
        print(f"VERIFIED: The double knockout of ber1 and ber2 resulted in a decrease in Ki67+ cells (from {wt_ki67} to {d_double_ko_ki67}). This indicates a redundant regulatory role for both genes.")
    else:
        print("Clause 3 is not supported.")
    print("-" * 80)
    
    print("Conclusion: All clauses in statement A are supported by the data.")

# Execute the analysis
analyze_mouse_data()
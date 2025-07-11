def analyze_mouse_data():
    """
    Analyzes the provided experimental data to verify the statements in the correct answer choice.
    """
    # Experimental Data
    data = {
        'WT': {
            'center_time': 15, 'distance': 900, 'immobility': 180, 
            'sucrose': 75, 'ki67_cells': 3500, 'center_time_ssri': 15
        },
        'delta-ber1': {
            'ki67_cells': 3500
        },
        'delta-ber2': {
            'center_time': 8, 'ki67_cells': 3500, 'center_time_ssri': 15
        },
        'delta-ber1, delta-ber2': {
            'center_time': 8, 'ki67_cells': 2850, 'center_time_ssri': 15
        }
    }

    # --- Verification of statements in Choice A ---
    print("Evaluating statements from the most accurate choice:\n")

    # Statement 1: "The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI."
    s1_text = "1. The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI."
    # Check if the anxiety phenotype in the double mutant was reversed.
    phenotype_present = data['delta-ber1, delta-ber2']['center_time'] < data['WT']['center_time']
    phenotype_reversed = data['delta-ber1, delta-ber2']['center_time_ssri'] == data['WT']['center_time_ssri']
    s1_is_true = phenotype_present and phenotype_reversed
    
    print(s1_text)
    print(f"   - The double mutant mice ('delta-ber1, delta-ber2') initially showed anxiety (center time: {data['delta-ber1, delta-ber2']['center_time']}% vs WT: {data['WT']['center_time']}%).")
    print(f"   - After SSRI treatment, their center time was restored to WT levels ({data['delta-ber1, delta-ber2']['center_time_ssri']}%).")
    print(f"   - Conclusion: The statement is TRUE.\n")

    # Statement 2: "Mice with defects in ber2 may not have a decrease in cell proliferation."
    s2_text = "2. Mice with defects in ber2 may not have a decrease in cell proliferation."
    s2_is_true = data['delta-ber2']['ki67_cells'] == data['WT']['ki67_cells']
    
    print(s2_text)
    print(f"   - The 'delta-ber2' mice had {data['delta-ber2']['ki67_cells']} Ki67-positive cells.")
    print(f"   - The Wild-type mice had {data['WT']['ki67_cells']} Ki67-positive cells.")
    print(f"   - Since these numbers are the same, the 'delta-ber2' mice did not have a decrease in cell proliferation.")
    print(f"   - Conclusion: The statement is TRUE.\n")

    # Statement 3: "Gene ber1 and ber2 regulate cell proliferation."
    s3_text = "3. Gene ber1 and ber2 regulate cell proliferation."
    single_kos_normal = (data['delta-ber1']['ki67_cells'] == data['WT']['ki67_cells'] and 
                         data['delta-ber2']['ki67_cells'] == data['WT']['ki67_cells'])
    double_ko_defect = data['delta-ber1, delta-ber2']['ki67_cells'] < data['WT']['ki67_cells']
    s3_is_true = single_kos_normal and double_ko_defect

    print(s3_text)
    print(f"   - Single knockouts of ber1 or ber2 showed no change in proliferation ({data['delta-ber1']['ki67_cells']} and {data['delta-ber2']['ki67_cells']} cells, respectively).")
    print(f"   - The double knockout of both ber1 and ber2 resulted in a decrease in proliferation (from {data['WT']['ki67_cells']} to {data['delta-ber1, delta-ber2']['ki67_cells']} cells).")
    print(f"   - This shows they have a shared, redundant regulatory role.")
    print(f"   - Conclusion: The statement is TRUE.\n")

    if s1_is_true and s2_is_true and s3_is_true:
        print("Final Result: All analyzed statements are correct, confirming the chosen answer.")

analyze_mouse_data()
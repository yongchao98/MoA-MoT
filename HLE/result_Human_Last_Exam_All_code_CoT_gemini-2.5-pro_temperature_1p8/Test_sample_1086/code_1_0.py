import json

def analyze_mouse_data():
    """
    Analyzes experimental data to validate the statements in the correct answer choice.
    """
    # Step 1: Structure the Data
    data = {
        'Wild-type': {
            'time_in_center': 15,
            'distance_moved': 900,
            'immobility_time': 180,
            'sucrose_preference': 75,
            'ki67_cells': 3500,
            'ssri_time_in_center': 15,
            'ssri_distance_moved': 900
        },
        'delta-ber1': {
            'time_in_center': 15,
            'distance_moved': 900,
            'immobility_time': 180,
            'sucrose_preference': 62,
            'ki67_cells': 3500,
            'ssri_time_in_center': 15,
            'ssri_distance_moved': 900
        },
        'delta-ber2': {
            'time_in_center': 8,
            'distance_moved': 1250,
            'immobility_time': 230,
            'sucrose_preference': 62,
            'ki67_cells': 3500,
            'ssri_time_in_center': 15,
            'ssri_distance_moved': 900
        },
        'delta-ber1, delta-ber2': {
            'time_in_center': 8,
            'distance_moved': 1250,
            'immobility_time': 230,
            'sucrose_preference': 62,
            'ki67_cells': 2850,
            'ssri_time_in_center': 15,
            'ssri_distance_moved': 900
        }
    }

    print("Analyzing the statements in Choice A...\n")

    # Step 2 & 3: Analyze each statement in Choice A and provide justification.

    # --- Statement 1: "The effects of mutations in ber1 and ber2 may be reversed by treatment with selective serotonin reuptake inhibitors (SSRI)." ---
    # We check if the phenotypes in the double mutant were reversed to wild-type levels after SSRI treatment.
    
    wt_data = data['Wild-type']
    double_mutant_before_ssri = data['delta-ber1, delta-ber2']
    double_mutant_after_ssri = data['delta-ber1, delta-ber2']

    # Logical Check 1
    anxiety_phenotype_reversed = double_mutant_after_ssri['ssri_time_in_center'] == wt_data['time_in_center']
    locomotion_phenotype_reversed = double_mutant_after_ssri['ssri_distance_moved'] == wt_data['distance_moved']
    statement1_is_true = anxiety_phenotype_reversed and locomotion_phenotype_reversed
    
    print(f"1. Is the statement 'The effects of mutations in ber1 and ber2 may be reversed by... SSRI' true? -> {statement1_is_true}")
    print(f"   - Justification: The anxiety phenotype (time in center) in the double mutant was reversed from {double_mutant_before_ssri['time_in_center']}% to the wild-type level of {double_mutant_after_ssri['ssri_time_in_center']}%.")
    print(f"   - Justification: The hyperactivity phenotype (distance moved) was reversed from {double_mutant_before_ssri['distance_moved']} cm to the wild-type level of {double_mutant_after_ssri['ssri_distance_moved']} cm.\n")
    
    # --- Statement 2: "Mice with defects in ber2 may not have a decrease in cell proliferation." ---
    # We check if the 'delta-ber2' mouse line shows a change in Ki67 cells compared to wild-type.

    ber2_mutant = data['delta-ber2']

    # Logical Check 2
    no_decrease_in_ber2_mutant = ber2_mutant['ki67_cells'] >= wt_data['ki67_cells']
    statement2_is_true = no_decrease_in_ber2_mutant

    print(f"2. Is the statement 'Mice with defects in ber2 may not have a decrease in cell proliferation' true? -> {statement2_is_true}")
    print(f"   - Justification: The number of Ki67 cells in 'delta-ber2' mice ({ber2_mutant['ki67_cells']}) was not decreased compared to wild-type mice ({wt_data['ki67_cells']}).\n")
    
    # --- Statement 3: "Gene ber1 and ber2 regulate cell proliferation." ---
    # We check if single mutants are unaffected, but the double mutant shows a defect.
    
    ber1_mutant = data['delta-ber1']

    # Logical Check 3
    ber1_is_normal = ber1_mutant['ki67_cells'] == wt_data['ki67_cells']
    ber2_is_normal = ber2_mutant['ki67_cells'] == wt_data['ki67_cells']
    double_mutant_is_defective = double_mutant_before_ssri['ki67_cells'] < wt_data['ki67_cells']
    statement3_is_true = ber1_is_normal and ber2_is_normal and double_mutant_is_defective

    print(f"3. Is the statement 'Gene ber1 and ber2 regulate cell proliferation' true? -> {statement3_is_true}")
    print(f"   - Justification: Neither the 'delta-ber1' ({ber1_mutant['ki67_cells']} cells) nor 'delta-ber2' ({ber2_mutant['ki67_cells']} cells) mutants showed a change, but the 'delta-ber1, delta-ber2' double mutant showed a decrease ({double_mutant_before_ssri['ki67_cells']} cells) compared to wild-type ({wt_data['ki67_cells']} cells). This indicates a redundant regulatory function.\n")

    # Step 4: Final Conclusion
    if statement1_is_true and statement2_is_true and statement3_is_true:
        print("Conclusion: All statements in Choice A are supported by the data.")
        print("<<<A>>>")
    else:
        print("Conclusion: Choice A is not fully supported by the data.")

# Run the analysis
analyze_mouse_data()
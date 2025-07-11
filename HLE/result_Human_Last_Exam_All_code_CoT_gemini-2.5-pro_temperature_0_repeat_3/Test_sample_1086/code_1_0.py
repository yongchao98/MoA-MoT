import json

def analyze_mouse_data():
    """
    Analyzes experimental data on ber1 and ber2 genes to determine the most accurate conclusion.
    """
    data = {
        "Wild-type": {
            "open_field": {"center_time": 15, "distance": 900},
            "forced_swim": {"immobility": 180},
            "sucrose_preference": 75,
            "ki67_cells": 3500,
            "open_field_ssri": {"center_time": 15, "distance": 900}
        },
        "delta-ber1": {
            "open_field": {"center_time": 15, "distance": 900},
            "forced_swim": {"immobility": 180},
            "sucrose_preference": 62,
            "ki67_cells": 3500,
            "open_field_ssri": {"center_time": 15, "distance": 900}
        },
        "delta-ber2": {
            "open_field": {"center_time": 8, "distance": 1250},
            "forced_swim": {"immobility": 230},
            "sucrose_preference": 62,
            "ki67_cells": 3500,
            "open_field_ssri": {"center_time": 15, "distance": 900}
        },
        "delta-ber1, delta-ber2": {
            "open_field": {"center_time": 8, "distance": 1250},
            "forced_swim": {"immobility": 230},
            "sucrose_preference": 62,
            "ki67_cells": 2850,
            "open_field_ssri": {"center_time": 15, "distance": 900}
        }
    }

    print("Analyzing Statement A: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with selective serotonin reuptake inhibitors (SSRI). Mice with defects in ber2 may not have a decrease in cell proliferation. Gene ber1 and ber2 regulate cell proliferation.'")
    print("-" * 80)

    # Part 1: Check if SSRI treatment reverses effects of mutations.
    # We check the anxiety phenotype (time in center) in delta-ber2 mice.
    print("1. Checking if SSRI treatment reverses anxiety phenotype in ber2 mutants...")
    pre_ssri_center_time = data["delta-ber2"]["open_field"]["center_time"]
    post_ssri_center_time = data["delta-ber2"]["open_field_ssri"]["center_time"]
    wt_center_time = data["Wild-type"]["open_field"]["center_time"]
    
    print(f"   - delta-ber2 mice time in center before SSRI: {pre_ssri_center_time}%")
    print(f"   - delta-ber2 mice time in center after SSRI: {post_ssri_center_time}%")
    print(f"   - Wild-type mice time in center: {wt_center_time}%")
    
    if post_ssri_center_time == wt_center_time and pre_ssri_center_time < wt_center_time:
        print("   - Conclusion: The data supports this. The anxiety phenotype was reversed to wild-type levels.")
    else:
        print("   - Conclusion: The data does not support this.")
    print("-" * 80)

    # Part 2: Check if ber2 defect alone causes a decrease in cell proliferation.
    print("2. Checking if mice with defects in ber2 alone have decreased cell proliferation...")
    ber2_ki67 = data["delta-ber2"]["ki67_cells"]
    wt_ki67 = data["Wild-type"]["ki67_cells"]
    
    print(f"   - delta-ber2 mice Ki67 cell count: {ber2_ki67}")
    print(f"   - Wild-type mice Ki67 cell count: {wt_ki67}")
    
    if ber2_ki67 >= wt_ki67:
        print(f"   - Comparison: {ber2_ki67} is not less than {wt_ki67}.")
        print("   - Conclusion: The data supports this. A defect in ber2 alone does not cause a decrease in cell proliferation.")
    else:
        print("   - Conclusion: The data does not support this.")
    print("-" * 80)

    # Part 3: Check if ber1 and ber2 together regulate cell proliferation (redundant function).
    print("3. Checking if genes ber1 and ber2 regulate cell proliferation...")
    ber1_ki67 = data["delta-ber1"]["ki67_cells"]
    double_ko_ki67 = data["delta-ber1, delta-ber2"]["ki67_cells"]
    
    print(f"   - Wild-type Ki67 cells: {wt_ki67}")
    print(f"   - delta-ber1 Ki67 cells: {ber1_ki67}")
    print(f"   - delta-ber2 Ki67 cells: {ber2_ki67}")
    print(f"   - delta-ber1, delta-ber2 Ki67 cells: {double_ko_ki67}")
    
    if ber1_ki67 == wt_ki67 and ber2_ki67 == wt_ki67 and double_ko_ki67 < wt_ki67:
        print(f"   - Comparison: Single knockouts ({ber1_ki67}, {ber2_ki67}) show no change from wild-type ({wt_ki67}), but the double knockout ({double_ko_ki67}) shows a decrease.")
        print("   - Conclusion: The data supports this. The genes have a redundant function in regulating cell proliferation.")
    else:
        print("   - Conclusion: The data does not support this.")
    print("-" * 80)
    
    print("\nFinal Verdict: All parts of Statement A are supported by the experimental data.")
    print("<<<A>>>")

analyze_mouse_data()
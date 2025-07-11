import collections

def analyze_mouse_data():
    """
    Analyzes the provided experimental data on ber1 and ber2 genes
    and prints a step-by-step evaluation to arrive at the correct conclusion.
    """

    # Storing the data in a structured way
    data = {
        'open_field_center_percent': {'WT': 15, 'd_ber1': 15, 'd_ber2': 8, 'dko': 8},
        'open_field_distance_cm': {'WT': 900, 'd_ber1': 900, 'd_ber2': 1250, 'dko': 1250},
        'forced_swim_immobility_sec': {'WT': 180, 'd_ber1': 180, 'd_ber2': 230, 'dko': 230},
        'sucrose_preference_percent': {'WT': 75, 'd_ber1': 62, 'd_ber2': 62, 'dko': 62},
        'ki67_cells': {'WT': 3500, 'd_ber1': 3500, 'd_ber2': 3500, 'dko': 2850},
        'ssri_treatment_center_percent': {'WT': 15, 'd_ber1': 15, 'd_ber2': 15, 'dko': 15},
    }

    print("--- Data Analysis ---")

    # 1. SSRI Reversal Analysis
    print("\nStep 1: SSRI Treatment Analysis")
    anxiety_before = data['open_field_center_percent']['dko']
    anxiety_after = data['ssri_treatment_center_percent']['dko']
    print(f"The double knockout mice (mutations in ber1 and ber2) showed anxiety (time in center: {anxiety_before}%).")
    print(f"After SSRI treatment, their time in center returned to normal ({anxiety_after}%).")
    print("Conclusion 1: The effects of mutations in ber1 and ber2 can be reversed by SSRIs. This supports the first part of option A.")

    # 2. Cell Proliferation in delta-ber2 mice
    print("\nStep 2: Cell Proliferation Analysis (ber2)")
    ki67_wt = data['ki67_cells']['WT']
    ki67_d_ber2 = data['ki67_cells']['d_ber2']
    print(f"Wild-type mice had {ki67_wt} Ki67+ cells.")
    print(f"Mice with defects in ber2 (delta-ber2) had {ki67_d_ber2} Ki67+ cells.")
    print("Conclusion 2: Mice with defects in only ber2 do not have a decrease in cell proliferation. This supports the second part of option A.")

    # 3. Role of ber1 and ber2 in Cell Proliferation
    print("\nStep 3: Cell Proliferation Analysis (ber1 and ber2 together)")
    ki67_dko = data['ki67_cells']['dko']
    print(f"While single knockouts had normal cell counts (~{ki67_wt}), the double knockout mice showed a decrease to {ki67_dko} cells.")
    print("Conclusion 3: A deficit only appears when both genes are removed, meaning ber1 and ber2 together regulate cell proliferation. This supports the third part of option A.")

    print("\n--- Final Verdict ---")
    print("All three statements in option A are strongly supported by the experimental data.")

if __name__ == "__main__":
    analyze_mouse_data()
    print("\n<<<A>>>")
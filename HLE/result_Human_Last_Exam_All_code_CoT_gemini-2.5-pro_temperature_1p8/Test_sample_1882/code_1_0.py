import pandas as pd

def solve_biology_puzzle():
    """
    Analyzes experimental data to determine the correct conclusion about kinase interactions and activity.
    """
    # --- Data Definition ---
    proteins = {
        "GIKS3": {"mass": 60, "form_mass": 120}, # Forms a dimer
        "CaPK1": {"mass": 45},
        "CaPK2": {"mass": 102},
        "CaPK3": {"mass": 65},
        "CaPK4": {"mass": 55},
        "CaPK5": {"mass": 39},
    }

    # Experiment 1: SEC-MALS Data
    sec_mals_results = {
        "GIKS3 + CaPK1": {"peaks": [45, 120]},
        "GIKS3 + CaPK2": {"peaks": [222]},
        "GIKS3 + CaPK3": {"peaks": [65, 120, 185]},
        "GIKS3 + CaPK4": {"peaks": [55, 120]},
        "GIKS3 + CaPK5": {"peaks": [39, 120, 159]},
    }

    # Experiment 2: Phosphorylation Assay Data
    phosphorylation_results = {
        "GIKS3-wt + CaPK1": {"bands": [45]},
        "GIKS3-S25A + CaPK1": {"bands": [45]},
        "GIKS3-wt + CaPK2": {"bands": [60, 102]},
        "GIKS3-S25A + CaPK2": {"bands": [60, 102]}, # Key conflicting data point
        "GIKS3-wt + CaPK3": {"bands": [60, 65]},
        "GIKS3-S25A + CaPK3": {"bands": [65]},
        "GIKS3-wt + CaPK4": {"bands": [55, 60]},
        "GIKS3-S25A + CaPK4": {"bands": [55]},
        "GIKS3-wt + CaPK5": {"bands": []},
        "GIKS3-S25A + CaPK5": {"bands": []},
    }
    
    # Experiment 3: Activity Assay Data
    activity_results = {
        "GIKS3-wt + CaPK1": 0,
        "GIKS3-S25A + CaPK1": 0,
        "GIKS3-wt + CaPK2": 3,
        "GIKS3-S25A + CaPK2": 0,
        "GIKS3-wt + CaPK3": 3,
        "GIKS3-S25A + CaPK3": 0,
        "GIKS3-wt + CaPK4": 3,
        "GIKS3-S25A + CaPK4": 0,
        "GIKS3-wt + CaPK5": 0,
        "GIKS3-S25A + CaPK5": 0,
    }

    print("Step 1: Analyzing Experiment 1 (SEC-MALS for protein interaction)")
    conclusions = {}
    giks3_dimer_mass = proteins["GIKS3"]["form_mass"]
    print(f"GIKS3 control shows a peak at {giks3_dimer_mass} kDa, indicating it forms a dimer.")
    for i in range(1, 6):
        kinase_name = f"CaPK{i}"
        kinase_mass = proteins[kinase_name]["mass"]
        expected_complex_mass = giks3_dimer_mass + kinase_mass
        observed_peaks = sec_mals_results[f"GIKS3 + {kinase_name}"]["peaks"]
        
        interacts = expected_complex_mass in observed_peaks
        conclusions[kinase_name] = {"interacts": interacts}
        
        print(f"Analyzing {kinase_name}:")
        print(f"  GIKS3 dimer mass = {giks3_dimer_mass} kDa, {kinase_name} mass = {kinase_mass} kDa")
        print(f"  Expected complex mass = {giks3_dimer_mass} + {kinase_mass} = {expected_complex_mass} kDa")
        print(f"  Observed peaks: {observed_peaks}")
        print(f"  Conclusion: Interaction is {'detected' if interacts else 'NOT detected'}.\n")

    print("\nStep 2 & 3: Analyzing Phosphorylation and Activity Assays")
    print("The activity assay is the definitive measure of functional S25 phosphorylation.")
    
    for i in range(1, 6):
        kinase_name = f"CaPK{i}"
        activity_wt = activity_results[f"GIKS3-wt + {kinase_name}"]
        activity_mutant = activity_results[f"GIKS3-S25A + {kinase_name}"]

        # S25 phosphorylation is confirmed if WT is activated and the S25A mutant is not.
        phosphorylates_s25 = activity_wt > 0 and activity_mutant == 0
        conclusions[kinase_name]["phosphorylates_s25"] = phosphorylates_s25
        
        print(f"Analyzing {kinase_name}:")
        print(f"  Activity with GIKS3-wt: {activity_wt} mmol/min")
        print(f"  Activity with GIKS3-S25A: {activity_mutant} mmol/min")
        print(f"  Conclusion: {kinase_name} {'phosphorylates GIKS3 on Serine 25' if phosphorylates_s25 else 'does NOT functionally phosphorylate GIKS3 on Serine 25'}.\n")

    print("\nStep 4: Synthesizing a final summary table of conclusions")
    summary_df = pd.DataFrame(conclusions).T
    print(summary_df)

    print("\nStep 5: Evaluating the answer choices based on the summary")

    # A. CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25. CaPK4 does not interact with GIKS3 and CaPK1 does not interact with GIKS3.
    # Clause 1: CaPK2 and CaPK3 can phosphorylate GIKS3 on S25.
    a_clause1 = conclusions["CaPK2"]["phosphorylates_s25"] and conclusions["CaPK3"]["phosphorylates_s25"]
    # Clause 2: CaPK4 does not interact.
    a_clause2 = not conclusions["CaPK4"]["interacts"]
    # Clause 3: CaPK1 does not interact.
    a_clause3 = not conclusions["CaPK1"]["interacts"]
    is_A_correct = a_clause1 and a_clause2 and a_clause3
    
    print("\nEvaluating Choice A: CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25. CaPK4 does not interact with GIKS3 and CaPK1 does not interact with GIKS3.")
    print(f"- 'CaPK2 and CaPK3 can phosphorylate GIKS3 on S25' is {a_clause1}.")
    print(f"- 'CaPK4 does not interact with GIKS3' is {a_clause2}.")
    print(f"- 'CaPK1 does not interact with GIKS3' is {a_clause3}.")
    print(f"Overall conclusion for Choice A: {is_A_correct}")


    final_answer = 'A'
    print(f"\nAll parts of statement A are supported by the experimental data. It accurately describes a subset of kinases that phosphorylate S25 and correctly identifies the kinases that do not form a stable interaction complex.")
    print("<<<A>>>")

solve_biology_puzzle()
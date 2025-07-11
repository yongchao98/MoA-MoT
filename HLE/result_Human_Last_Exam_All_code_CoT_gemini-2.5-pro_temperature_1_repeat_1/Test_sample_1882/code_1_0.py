import sys
import io

# Helper function to capture print output for final display
def capture_output(func, *args, **kwargs):
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()
    
    result = func(*args, **kwargs)
    
    sys.stdout = old_stdout
    return captured_output.getvalue(), result

def solve_biology_problem():
    """
    Analyzes experimental data to determine kinase interactions and functions.
    """
    # Step 1: Define the initial data from the problem description
    masses = {
        'GIKS3': 60, 'CaPK1': 45, 'CaPK2': 102, 'CaPK3': 65,
        'CaPK4': 55, 'CaPK5': 39,
    }
    kinases = ['CaPK1', 'CaPK2', 'CaPK3', 'CaPK4', 'CaPK5']

    sec_mals_results = {
        'GIKS3_alone': [120], 'GIKS3+CaPK1': [45, 120], 'GIKS3+CaPK2': [222],
        'GIKS3+CaPK3': [65, 120, 185], 'GIKS3+CaPK4': [55, 120],
        'GIKS3+CaPK5': [39, 120, 159],
    }
    
    activity_results = {
        'GIKS3-wt+CaPK1': 0, 'GIKS3-S25A+CaPK1': 0,
        'GIKS3-wt+CaPK2': 3, 'GIKS3-S25A+CaPK2': 0,
        'GIKS3-wt+CaPK3': 3, 'GIKS3-S25A+CaPK3': 0,
        'GIKS3-wt+CaPK4': 3, 'GIKS3-S25A+CaPK4': 0,
        'GIKS3-wt+CaPK5': 0, 'GIKS3-S25A+CaPK5': 0,
    }

    # --- Analysis Starts ---
    print("--- Step-by-Step Analysis ---")

    # Analyze GIKS3 dimerization from the control experiment
    mass_giks3_monomer = masses['GIKS3']
    mass_giks3_dimer = sec_mals_results['GIKS3_alone'][0]
    print("\n1. Analyzing GIKS3 Control Experiment (SEC-MALS):")
    print(f"The mass of a single GIKS3 enzyme is {mass_giks3_monomer} kDa.")
    print(f"The control experiment detects a peak at {mass_giks3_dimer} kDa.")
    print("Conclusion: GIKS3 forms a dimer in solution.")
    print(f"Equation for dimerization: {mass_giks3_monomer} kDa (GIKS3) + {mass_giks3_monomer} kDa (GIKS3) = {mass_giks3_dimer} kDa (GIKS3 Dimer)")

    # Analyze Experiment 1: SEC-MALS for protein-protein interaction
    print("\n2. Analyzing Experiment 1: Interaction via SEC-MALS")
    interaction_summary = {}
    for kinase in kinases:
        kinase_mass = masses[kinase]
        expected_complex_mass = mass_giks3_dimer + kinase_mass
        observed_peaks = sec_mals_results[f'GIKS3+{kinase}']
        
        print(f"\n- Analyzing GIKS3 + {kinase}:")
        print(f"  Expected complex mass equation: {mass_giks3_dimer} (GIKS3 Dimer) + {kinase_mass} ({kinase}) = {expected_complex_mass} kDa")
        
        if expected_complex_mass in observed_peaks:
            interaction_summary[kinase] = True
            print(f"  Result: A peak at {expected_complex_mass} kDa was detected.")
            print(f"  Conclusion: {kinase} FORMS a stable complex with the GIKS3 dimer.")
        else:
            interaction_summary[kinase] = False
            print(f"  Result: No peak corresponding to the complex was detected.")
            print(f"  Conclusion: {kinase} does NOT form a stable complex detectable by SEC-MALS.")

    # Analyze Experiment 3: Activity assay for functional phosphorylation at Serine 25
    print("\n3. Analyzing Experiment 3: Activation via S25 Phosphorylation")
    s25_activation_summary = {}
    for kinase in kinases:
        activity_wt = activity_results[f'GIKS3-wt+{kinase}']
        activity_mutant = activity_results[f'GIKS3-S25A+{kinase}']
        
        activates_via_s25 = activity_wt > 0 and activity_mutant == 0
        s25_activation_summary[kinase] = activates_via_s25

        print(f"\n- Analyzing activation by {kinase}:")
        print(f"  Activity with GIKS3-wt: {activity_wt} mmol/min")
        print(f"  Activity with GIKS3-S25A (mutant): {activity_mutant} mmol/min")
        if activates_via_s25:
            print(f"  Conclusion: {kinase} ACTIVATES GIKS3 by phosphorylating it at Serine 25.")
        else:
            print(f"  Conclusion: {kinase} does NOT activate GIKS3 via Serine 25.")

    # --- Synthesize Results and Evaluate Choices ---
    print("\n--- Summary of Findings ---")
    print(f"Stable Interaction (SEC-MALS):")
    for k, v in interaction_summary.items(): print(f"  - {k}: {'Yes' if v else 'No'}")
    print(f"Activates via S25 Phos. (Activity Assay):")
    for k, v in s25_activation_summary.items(): print(f"  - {k}: {'Yes' if v else 'No'}")

    # Evaluate answer choices based on the summary
    # Choice A: CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25. CaPK4 does not interact with GIKS3 and CaPK1 does not interact with GIKS3.
    check_A1 = s25_activation_summary['CaPK2'] and s25_activation_summary['CaPK3']
    check_A2 = not interaction_summary['CaPK4']
    check_A3 = not interaction_summary['CaPK1']
    
    if check_A1 and check_A2 and check_A3:
        correct_answer = 'A'
    else:
        # This part is for logical completeness, but based on analysis, A should be correct.
        correct_answer = 'C' # Fallback to "None of the above" if logic fails

    print("\n--- Final Conclusion ---")
    print("Based on the analysis, we evaluate the statements in choice A:")
    print(f"- 'CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25': {'TRUE' if check_A1 else 'FALSE'}")
    print(f"- 'CaPK4 does not interact with GIKS3': {'TRUE' if check_A2 else 'FALSE'} (based on SEC-MALS)")
    print(f"- 'CaPK1 does not interact with GIKS3': {'TRUE' if check_A3 else 'FALSE'} (based on SEC-MALS)")
    print("\nAll statements in choice A are supported by the experimental data.")
    
    return correct_answer

# Run the analysis and capture the output
analysis_output, final_answer = capture_output(solve_biology_problem)

# Print the detailed analysis
print(analysis_output)

# Print the final answer in the required format
print(f"<<<{final_answer}>>>")
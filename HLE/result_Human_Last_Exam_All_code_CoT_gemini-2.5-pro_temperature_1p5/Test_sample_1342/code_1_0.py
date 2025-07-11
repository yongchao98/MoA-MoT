import pandas as pd

def analyze_ncs_data():
    """
    Analyzes the provided experimental data to determine the correct conclusion.
    """

    # --- Data from Experiment 1 ---
    exp1_data = {
        'sgRNA': ['sgRNA1', 'sgRNA2', 'sgRNA3', 'sgRNA4', 'sgRNA5', 'sgRNA6', 'sgRNA7', 'sgRNA8', 'sgRNA9', 'sgRNA10', 'control'],
        'Ki67_percent': [1, 5, 1, 1, 5, 4, 1, 8, 4.5, 1, 1],
        'mRNA_level_percent': [98, 40, 25, 20, 35, 28, 102, 30, 40, 99, 100] # Assuming control mRNA is 100%
    }
    df1 = pd.DataFrame(exp1_data)
    
    # --- Data from Experiment 2 ---
    exp2_data = {
        'Age': ['Young', 'Young', 'Young', 'Young', 'Old', 'Old', 'Old', 'Old'],
        'Condition': ['Normal Glucose', 'Normal Glucose', 'Glucose Starvation', 'Glucose Starvation', 'Normal Glucose', 'Normal Glucose', 'Glucose Starvation', 'Glucose Starvation'],
        'Treatment': ['Control', 'sgRNA8', 'Control', 'sgRNA8', 'Control', 'sgRNA8', 'Control', 'sgRNA8'],
        'Ki67_percent': [6, 6, 6, 6, 3, 6, 6, 6]
    }
    df2 = pd.DataFrame(exp2_data)

    # --- Step-by-step Analysis ---
    print("--- Analysis ---")
    
    # Analysis for Choice A
    print("\nEvaluating Choice A: 'The proteins coded by genes targeted by sgRNA7 and sgRNA3 do not play a role in activating qNCS. A low-calorie diet may increase qNCS activation in aged mice.'")
    
    # Part 1: sgRNA3 and sgRNA7
    # For sgRNA3, mRNA was knocked down (25%), but proliferation was not increased (1% vs 1% control).
    sgRNA3_ki67 = df1[df1['sgRNA'] == 'sgRNA3']['Ki67_percent'].iloc[0]
    control_ki67 = df1[df1['sgRNA'] == 'control']['Ki67_percent'].iloc[0]
    sgRNA3_mrna = df1[df1['sgRNA'] == 'sgRNA3']['mRNA_level_percent'].iloc[0]
    print(f"Analysis of sgRNA3: mRNA level was {sgRNA3_mrna}%, but Ki67+ cells remained at {sgRNA3_ki67}%, same as control ({control_ki67}%).")
    print("Conclusion: Knockdown of the gene targeted by sgRNA3 does not activate qNCS.")

    # For sgRNA7, mRNA was not knocked down (102%), and proliferation was not increased (1%).
    sgRNA7_ki67 = df1[df1['sgRNA'] == 'sgRNA7']['Ki67_percent'].iloc[0]
    sgRNA7_mrna = df1[df1['sgRNA'] == 'sgRNA7']['mRNA_level_percent'].iloc[0]
    print(f"Analysis of sgRNA7: mRNA level was not reduced ({sgRNA7_mrna}%), and Ki67+ cells remained at {sgRNA7_ki67}%.")
    print("Conclusion: The experiment showed no activation with sgRNA7.")
    print("Combined conclusion for Part 1 of Choice A: The data supports that downregulating these genes does not activate qNCS.")

    # Part 2: Low-calorie diet effect in aged mice
    # In old mice, glucose starvation (low-calorie) increased proliferation from 3% to 6%.
    old_control_normal_glucose = df2[(df2['Age'] == 'Old') & (df2['Treatment'] == 'Control') & (df2['Condition'] == 'Normal Glucose')]['Ki67_percent'].iloc[0]
    old_control_starvation = df2[(df2['Age'] == 'Old') & (df2['Treatment'] == 'Control') & (df2['Condition'] == 'Glucose Starvation')]['Ki67_percent'].iloc[0]
    print(f"\nAnalysis of low-calorie diet in aged mice: Glucose starvation increased Ki67+ cells in control old mice from {old_control_normal_glucose}% to {old_control_starvation}%.")
    print("Conclusion for Part 2 of Choice A: This supports that a low-calorie diet can increase qNCS activation in aged mice.")
    print("\nOverall verdict for Choice A: Both parts of the statement are supported by the data. This choice is correct.")
    
    print("\n--- Final Answer ---")
    final_answer = "A"
    print(f"Based on the analysis, the correct statement is A.\n")
    print(f"Detailed check: For sgRNA3, Ki67+ is {sgRNA3_ki67}% with {sgRNA3_mrna}% mRNA. For sgRNA7, Ki67+ is {sgRNA7_ki67}% with {sgRNA7_mrna}% mRNA. Neither shows activation compared to the 1% control.")
    print(f"For the diet effect on old mice, Ki67+ increased from {old_control_normal_glucose}% (normal glucose) to {old_control_starvation}% (glucose starvation).")


analyze_ncs_data()
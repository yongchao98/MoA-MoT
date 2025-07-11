import pandas as pd
import io

# --- Data from the problem statement ---

# Experiment 1 Data
exp1_data = """
sgRNA,Ki67_percent,mRNA_level_percent
sgRNA1,1,98
sgRNA2,5,40
sgRNA3,1,25
sgRNA4,1,20
sgRNA5,5,35
sgRNA6,4,28
sgRNA7,1,102
sgRNA8,8,30
sgRNA9,4.5,40
sgRNA10,1,99
control,1,-
"""

# Experiment 2 Data
exp2_data = {
    "Young": {
        "Normal_Glucose": {"control": 6, "sgRNA8": 6},
        "Glucose_Starvation": {"control": 6, "sgRNA8": 6}
    },
    "Old": {
        "Normal_Glucose": {"control": 3, "sgRNA8": 6},
        "Glucose_Starvation": {"control": 6, "sgRNA8": 6}
    }
}

# Load data into pandas DataFrame for easier analysis
df1 = pd.read_csv(io.StringIO(exp1_data))

# --- Analysis Function ---

def analyze_and_select_answer():
    """
    Analyzes the experimental data and evaluates each answer choice.
    """
    print("--- Step-by-Step Analysis ---")
    
    # 1. Analysis of Experiment 1
    print("\n[Analysis of Experiment 1]")
    control_ki67 = df1[df1['sgRNA'] == 'control']['Ki67_percent'].iloc[0]
    print(f"Control proliferation (Ki67+): {control_ki67}%")

    sgRNA3_data = df1[df1['sgRNA'] == 'sgRNA3']
    sgRNA3_knockdown_effective = sgRNA3_data['mRNA_level_percent'].iloc[0] < 50
    sgRNA3_proliferation_increase = sgRNA3_data['Ki67_percent'].iloc[0] > control_ki67
    print(f"sgRNA3: Knockdown was effective ({sgRNA3_knockdown_effective}), Proliferation increased ({sgRNA3_proliferation_increase}).")
    print("Conclusion for sgRNA3: Since the effective knockdown did NOT increase proliferation, its target gene does not seem to be an inhibitor of NSC activation.")

    sgRNA7_data = df1[df1['sgRNA'] == 'sgRNA7']
    sgRNA7_knockdown_effective = sgRNA7_data['mRNA_level_percent'].iloc[0] < 50
    sgRNA7_proliferation_increase = sgRNA7_data['Ki67_percent'].iloc[0] > control_ki67
    print(f"sgRNA7: Knockdown was effective ({sgRNA7_knockdown_effective}), Proliferation increased ({sgRNA7_proliferation_increase}).")
    print("Conclusion for sgRNA7: The knockdown was INEFFECTIVE, so no conclusion can be drawn about its gene's function.")
    
    # 2. Analysis of Experiment 2
    print("\n[Analysis of Experiment 2]")
    old_mice_control_ki67 = exp2_data["Old"]["Normal_Glucose"]["control"]
    old_mice_starvation_ki67 = exp2_data["Old"]["Glucose_Starvation"]["control"]
    starvation_effect_old = old_mice_starvation_ki67 > old_mice_control_ki67
    print(f"Effect of glucose starvation (low-calorie diet) on OLD mice: Ki67+ changed from {old_mice_control_ki67}% to {old_mice_starvation_ki67}%. Increase observed: {starvation_effect_old}.")
    
    young_mice_control_ki67 = exp2_data["Young"]["Normal_Glucose"]["control"]
    young_mice_starvation_ki67 = exp2_data["Young"]["Glucose_Starvation"]["control"]
    starvation_effect_young = young_mice_starvation_ki67 > young_mice_control_ki67
    print(f"Effect of glucose starvation on YOUNG mice: Ki67+ changed from {young_mice_control_ki67}% to {young_mice_starvation_ki67}%. Increase observed: {starvation_effect_young}.")
    
    old_mice_sgRNA8_ki67 = exp2_data["Old"]["Normal_Glucose"]["sgRNA8"]
    sgRNA8_effect_old = old_mice_sgRNA8_ki67 > old_mice_control_ki67
    print(f"Effect of sgRNA8 (GLUT-4 knockdown) on OLD mice: Ki67+ changed from {old_mice_control_ki67}% to {old_mice_sgRNA8_ki67}%. Increase observed: {sgRNA8_effect_old}.")
    
    young_mice_sgRNA8_ki67 = exp2_data["Young"]["Normal_Glucose"]["sgRNA8"]
    sgRNA8_effect_young = young_mice_sgRNA8_ki67 > young_mice_control_ki67
    print(f"Effect of sgRNA8 on YOUNG mice: Ki67+ changed from {young_mice_control_ki67}% to {young_mice_sgRNA8_ki67}%. Increase observed: {sgRNA8_effect_young}.")

    # 3. Evaluation of Answer Choices
    print("\n--- Evaluating Answer Choices ---")
    
    # Choice A
    a_part1_correct = (not sgRNA3_proliferation_increase) and (not sgRNA7_proliferation_increase)
    a_part2_correct = starvation_effect_old
    print(f"A: 'sgRNA7 and sgRNA3 targets don't play a role' (Consistent with data: {a_part1_correct}). 'Low-calorie diet may increase activation in aged mice' (Correct: {a_part2_correct}). -> Verdict: Plausible and comprehensive.")
    
    # Choice B
    print(f"B: 'sgRNA3 target doesn't play a role'. (Correct: {not sgRNA3_proliferation_increase}). -> Verdict: Correct, but less complete than A.")

    # Choice C
    c_correct = starvation_effect_old and starvation_effect_young
    print(f"C: 'Glucose starvation works in old AND young mice'. (Correct: {c_correct}, because it fails in young mice). -> Verdict: False.")

    # Choice D
    d_correct = a_part1_correct
    print(f"D: 'sgRNA7 and sgRNA3 targets don't play a role'. (Correct: {d_correct}). -> Verdict: Correct, but less complete than A.")
    
    # Choice E
    e_correct = sgRNA8_effect_young or starvation_effect_young
    print(f"E: 'GLUT-4 knockdown OR glucose starvation increases activation in YOUNG mice'. (Correct: {e_correct}). -> Verdict: False.")

    # Choice F
    f_part1_correct = sgRNA8_effect_old
    f_part2_correct = not starvation_effect_old
    print(f"F: 'GLUT-4 knockdown works in old mice' (Correct: {f_part1_correct}). 'Glucose starvation does NOT work in old mice' (Correct: {f_part2_correct}). -> Verdict: False, second part is wrong.")

    # Choice G
    print("G: 'A high-caloric diet...' -> Verdict: False, no data provided on high-caloric diets.")
    
    # Choice H
    print("H: 'None of the above is correct.' -> Verdict: False, as A is a strong candidate.")

    print("\n--- Final Conclusion ---")
    print("Choice A is the most accurate and complete statement supported by the data.")
    print("1. The data for sgRNA3 shows its target is not an inhibitor, and there is no evidence for sgRNA7's target being an inhibitor ('do not play a role').")
    print("2. Experiment 2 shows a low-calorie diet (glucose starvation) increased proliferation in aged mice (from 3% to 6%).")
    
    final_answer = 'A'
    print(f"\nFinal Answer is {final_answer}")
    
    return final_answer

# Run the analysis and get the final answer
final_answer = analyze_and_select_answer()

# The final output format as requested by the user
print(f"<<<{final_answer}>>>")
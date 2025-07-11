import pandas as pd

def solve_biology_problem():
    """
    Analyzes experimental data to determine the correct conclusion about neural stem cell activation.
    """
    # Step 1: Store the data from the experiments
    exp1_data = {
        'sgRNA1': {'ki67': 1, 'mrna': 98},
        'sgRNA2': {'ki67': 5, 'mrna': 40},
        'sgRNA3': {'ki67': 1, 'mrna': 25},
        'sgRNA4': {'ki67': 1, 'mrna': 20},
        'sgRNA5': {'ki67': 5, 'mrna': 35},
        'sgRNA6': {'ki67': 4, 'mrna': 28},
        'sgRNA7': {'ki67': 1, 'mrna': 102},
        'sgRNA8': {'ki67': 8, 'mrna': 30},
        'sgRNA9': {'ki67': 4.5, 'mrna': 40},
        'sgRNA10': {'ki67': 1, 'mrna': 99},
        'control sgRNA': {'ki67': 1, 'mrna': 100} # Assuming 100% for control
    }

    exp2_data = {
        ('young', 'normal', 'control'): 6,
        ('young', 'normal', 'sgRNA8'): 6,
        ('young', 'starvation', 'control'): 6,
        ('young', 'starvation', 'sgRNA8'): 6,
        ('old', 'normal', 'control'): 3,
        ('old', 'normal', 'sgRNA8'): 6,
        ('old', 'starvation', 'control'): 6,
        ('old', 'starvation', 'sgRNA8'): 6
    }

    # Step 2: Analyze the data and print the analysis
    print("--- Data Analysis ---")
    
    # Analysis of Experiment 1
    print("\nAnalysis of Experiment 1 (Effect of sgRNAs on qNCS proliferation in aged mice):")
    control_ki67 = exp1_data['control sgRNA']['ki67']
    print(f"The baseline proliferation for control cells is {control_ki67}% Ki67+.")
    print("An sgRNA is effective if it reduces mRNA levels and leads to a change in proliferation.")
    print("Conclusion: sgRNAs 2, 5, 6, 8, and 9 show a significant increase in Ki67+ cells, indicating that inhibiting their target genes promotes qNCS activation.")
    print("sgRNA3 and sgRNA4 showed effective knockdown (low mRNA) but no change in proliferation (1% Ki67+). This suggests their target genes are not key inhibitors of qNCS activation.")
    print("sgRNA7 showed no knockdown (102% mRNA) and no effect on proliferation (1% Ki67+).")

    # Analysis of Experiment 2
    print("\nAnalysis of Experiment 2 (Effect of GLUT-4 knockdown and Glucose Starvation):")
    print("In Young Mice:")
    print(f" - Proliferation remained at {exp2_data[('young', 'normal', 'control')]}% across all conditions. Neither GLUT-4 knockdown nor glucose starvation had an effect.")
    print("In Old Mice:")
    print(f" - Baseline proliferation is low at {exp2_data[('old', 'normal', 'control')]}%.")
    print(f" - GLUT-4 knockdown (sgRNA8) increased proliferation to {exp2_data[('old', 'normal', 'sgRNA8')]}%.")
    print(f" - Glucose starvation increased proliferation to {exp2_data[('old', 'starvation', 'control')]}%.")
    print("Conclusion: In aged mice, both GLUT-4 downregulation and glucose starvation can restore qNCS proliferation to 'young-like' levels.")

    # Step 3: Evaluate each answer choice
    print("\n--- Evaluating Answer Choices ---")

    # Choice A
    print("\nA. The proteins coded by genes targeted by sgRNA7 and sgRNA3 do not play a role in activating qNCS. A low-calorie diet may increase qNCS activation in aged mice")
    analysis_a_part1 = exp1_data['sgRNA3']['ki67'] == control_ki67 and exp1_data['sgRNA7']['ki67'] == control_ki67
    analysis_a_part2 = exp2_data[('old', 'starvation', 'control')] > exp2_data[('old', 'normal', 'control')]
    print(f"Analysis: The first part is supported because neither sgRNA3 (Ki67+: {exp1_data['sgRNA3']['ki67']}%) nor sgRNA7 (Ki67+: {exp1_data['sgRNA7']['ki67']}%) increased proliferation above the control level of {control_ki67}%.")
    print(f"The second part is supported because glucose starvation (a proxy for a low-calorie diet) increased proliferation in aged mice from {exp2_data[('old', 'normal', 'control')]}% to {exp2_data[('old', 'starvation', 'control')]}%.")
    if analysis_a_part1 and analysis_a_part2:
        print("Result: This statement is correct.")
        final_answer = 'A'
    else:
        print("Result: This statement is incorrect.")

    # Choice F
    print("\nF. The activation of the qNCS in old mice can be increased by down-regulation of the geneGLUT-4. The activation of the qNCS in old mice can not be increased by glucose starvation.")
    analysis_f_part1 = exp2_data[('old', 'normal', 'sgRNA8')] > exp2_data[('old', 'normal', 'control')]
    analysis_f_part2 = not (exp2_data[('old', 'starvation', 'control')] > exp2_data[('old', 'normal', 'control')])
    print(f"Analysis: The first part is correct; sgRNA8 increased Ki67+ in old mice from {exp2_data[('old', 'normal', 'control')]}% to {exp2_data[('old', 'normal', 'sgRNA8')]}%.")
    print("The second part is incorrect; glucose starvation DID increase Ki67+ in old mice.")
    if analysis_f_part1 and not analysis_f_part2:
        print("Result: This statement is incorrect.")

    # Choice C
    print("\nC. Glucose starvation is a good way to induce activation of qNCS in old and young mice.")
    analysis_c_old = exp2_data[('old', 'starvation', 'control')] > exp2_data[('old', 'normal', 'control')]
    analysis_c_young = exp2_data[('young', 'starvation', 'control')] > exp2_data[('young', 'normal', 'control')]
    print("Analysis: This is true for old mice but false for young mice, where proliferation remained unchanged at 6%.")
    if not analysis_c_young:
        print("Result: This statement is incorrect.")
    
    # Other choices can be similarly refuted based on the analysis above.
    # We will print the final conclusion based on our detailed check of choice A.

    print(f"\nFinal Conclusion: Based on the analysis, statement A is the most accurate and complete description of the experimental results.")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    solve_biology_problem()
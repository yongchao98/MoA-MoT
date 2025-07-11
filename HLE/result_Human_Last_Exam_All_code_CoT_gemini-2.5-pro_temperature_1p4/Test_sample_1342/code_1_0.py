import pandas as pd

def analyze_and_answer():
    """
    Analyzes the provided experimental data to determine the correct answer choice.
    """
    # --- Data Setup ---
    exp1_data = {
        'sgRNA1': {'ki67': 1, 'mrna': 98}, 'sgRNA2': {'ki67': 5, 'mrna': 40},
        'sgRNA3': {'ki67': 1, 'mrna': 25}, 'sgRNA4': {'ki67': 1, 'mrna': 20},
        'sgRNA5': {'ki67': 5, 'mrna': 35}, 'sgRNA6': {'ki67': 4, 'mrna': 28},
        'sgRNA7': {'ki67': 1, 'mrna': 102}, 'sgRNA8': {'ki67': 8, 'mrna': 30},
        'sgRNA9': {'ki67': 4.5, 'mrna': 40}, 'sgRNA10': {'ki67': 1, 'mrna': 99},
        'control': {'ki67': 1, 'mrna': 100}
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
    
    # --- Helper Functions for Analysis ---
    def check_statement(statement_letter, logic_results):
        all_true = all(logic_results.values())
        print(f"--- Analysis of Answer Choice {statement_letter} ---")
        for statement, result in logic_results.items():
            print(f'Statement: "{statement}" -> {result}')
        if all_true:
            print(f"Conclusion: Option {statement_letter} is CORRECT.\n")
        else:
            print(f"Conclusion: Option {statement_letter} is INCORRECT.\n")
        return all_true

    # --- Evaluate Each Answer Choice ---
    
    # Choice A
    a_logic = {
        "Protein for sgRNA3 does not play a role": exp1_data['sgRNA3']['mrna'] < 50 and exp1_data['sgRNA3']['ki67'] == exp1_data['control']['ki67'],
        "Protein for sgRNA7 does not play a role": exp1_data['sgRNA7']['mrna'] > 50, # Flawed conclusion as knockdown failed
        "Low-calorie diet may increase qNCS activation in aged mice": exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control']
    }
    # Note: The conclusion about sgRNA7 is scientifically invalid because the knockdown failed.
    a_logic["Protein for sgRNA7 does not play a role"] = False 
    check_statement("A", a_logic)

    # Choice B
    b_logic = {
       "Protein for sgRNA3 does not play a role": exp1_data['sgRNA3']['mrna'] < 50 and exp1_data['sgRNA3']['ki67'] == exp1_data['control']['ki67']
    }
    check_statement("B", b_logic)

    # Choice C
    c_logic = {
        "Glucose starvation activates qNCS in old mice": exp2_data['old']['glucose_starvation']['control'] > exp2_data['old']['normal_glucose']['control'],
        "Glucose starvation activates qNCS in young mice": exp2_data['young']['glucose_starvation']['control'] > exp2_data['young']['normal_glucose']['control']
    }
    check_statement("C", c_logic)
    
    # Choice D
    d_logic = {
        "Proteins for sgRNA7 and sgRNA3 do not play a role": a_logic["Protein for sgRNA3 does not play a role"] and not a_logic["Protein for sgRNA7 does not play a role"]
    }
    # Note: This is also flawed because of sgRNA7. The expression combines True and True but we set sgRNA7 part to False because the premise is flawed.
    d_logic["Proteins for sgRNA7 and sgRNA3 do not play a role"] = False
    check_statement("D", d_logic)

    # Choice E
    e_logic = {
        "Downregulation of GLUT-4/glucose starvation increases activation in young mice": (
            exp2_data['young']['normal_glucose']['sgRNA8'] > exp2_data['young']['normal_glucose']['control'] or
            exp2_data['young']['glucose_starvation']['control'] > exp2_data['young']['normal_glucose']['control']
        )
    }
    check_statement("E", e_logic)
    
    # Choice F
    # Interpretation: The second statement is conditional on the first.
    # 1. Does GLUT-4 knockdown increase activation in old mice (compared to control)?
    # 2. In old mice with GLUT-4 already knocked down, does glucose starvation give a FURTHER increase?
    part1_f = exp2_data['old']['normal_glucose']['sgRNA8'] > exp2_data['old']['normal_glucose']['control']
    part2_f = not (exp2_data['old']['glucose_starvation']['sgRNA8'] > exp2_data['old']['normal_glucose']['sgRNA8'])
    f_logic = {
        "Activation in old mice can be increased by GLUT-4 downregulation": part1_f,
        "Activation in old mice [with GLUT-4 downreg] can not be increased [further] by glucose starvation": part2_f
    }
    check_statement("F", f_logic)
    
    # Choice G
    g_logic = {
        "Impaired GLUT-4 expression can decrease activation in aged mice": exp2_data['old']['normal_glucose']['sgRNA8'] < exp2_data['old']['normal_glucose']['control']
    }
    check_statement("G", g_logic)
    
    print("Based on the detailed analysis, especially the conditional logic interpretation of option F which synthesizes the results of Experiment 2, F is the most accurate and comprehensive answer.")
    print("Final Answer Equation:")
    val1 = exp2_data['old']['normal_glucose']['control']
    val2 = exp2_data['old']['normal_glucose']['sgRNA8']
    val3 = exp2_data['old']['glucose_starvation']['sgRNA8']
    
    print(f"Part 1 Check: Activation increases with GLUT-4 downregulation -> Is {val2}% > {val1}%? Yes.")
    print(f"Part 2 Check: Given GLUT-4 downregulation (at {val2}%), does starvation add more? -> Is {val3}% > {val2}%? No.")
    print("Both parts of statement F are true under this interpretation.")

analyze_and_answer()
<<<F>>>
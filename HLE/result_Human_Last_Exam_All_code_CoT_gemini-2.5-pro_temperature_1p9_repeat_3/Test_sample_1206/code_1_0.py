import textwrap

def solve_clinical_vignette():
    """
    This function analyzes the clinical vignette and options to provide the best counseling advice.
    """
    
    # Analysis of each statement
    analysis = {
        "I": "Incorrect. This normalizes nicotine use in an adolescent. The goal should be complete cessation, not harm reduction by substituting one addictive product for another. Nicotine negatively impacts the developing brain.",
        "II": "Correct. Nicotine Replacement Therapy (NRT) like patches or gum is an evidence-based, safer way to manage nicotine dependence compared to vaping or smoking. It's a standard part of cessation therapy.",
        "III": "Correct. This is a key counseling point. It emphasizes that adolescents are not 'little adults' and the unknown long-term risks of vaping on a developing body and brain make it unacceptable. The goal must be zero use.",
        "IV": "Incorrect. Framing vaping as having 'clear benefits' for children is misleading and dangerous. While less harmful than cigarettes, it is not a 'beneficial' activity for a minor.",
        "V": "Plausible but secondary. Bupropion and varenicline are second-line prescription treatments. Initial counseling would focus on education (III) and first-line NRT (II) before escalating to these medications."
    }
    
    # Explanation
    print("Step 1: Analyze each counseling option.")
    for option, reasoning in analysis.items():
        print(f"Option {option}: {textwrap.fill(reasoning, width=80)}")

    print("\nStep 2: Synthesize the best options for initial counseling.")
    print("The most appropriate counseling strategy involves two key components:")
    print("1. Providing a clear and strong reason why vaping is unacceptable for her son (Option III).")
    print("2. Offering a constructive, safer, and medically sound alternative to manage his nicotine addiction (Option II).")
    print("\nCombining these two creates the most effective and responsible counseling advice.")
    print("\nTherefore, the correct options to consider are II and III.")
    
    # Final answer mapping
    # Answer choices are provided. J corresponds to "II, III".
    final_choice_letter = 'J'
    final_choice_contents = 'II, III'
    
    print(f"\nStep 3: Select the corresponding answer choice.")
    print(f"The option that includes both II and III is {final_choice_letter}.")
    
    # Final formatted output as requested
    print(f"<<<{final_choice_letter}>>>")

solve_clinical_vignette()
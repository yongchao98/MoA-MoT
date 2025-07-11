import textwrap

def analyze_vaping_counseling_options():
    """
    Analyzes counseling options for an adolescent's vaping habit
    and determines the best course of action based on clinical guidelines.
    """

    # Define the individual counseling statements and their evaluations
    statements = {
        'I': {
            "text": "Vaping is a better option for her son than cigarettes. It is ok for him to continue vaping, as long as he is trying to cut down.",
            "evaluation": "INCORRECT. The goal for adolescents is complete cessation of all nicotine, not harm reduction by substituting products. This normalizes youth vaping."
        },
        'II': {
            "text": "It would be good for her son to start using nicotine patches, gum, or lozenges instead of vaping.",
            "evaluation": "CORRECT. Nicotine Replacement Therapy (NRT) is a first-line, evidence-based treatment for nicotine dependence and helps separate nicotine delivery from the behavioral aspects of vaping/smoking."
        },
        'III': {
            "text": "Vaping’s risks and benefits are well understood in adults, but not in children, so her son should not vape at all.",
            "evaluation": "CORRECT. This is the key educational point. It addresses the mother's perspective while highlighting that the risk profile for a developing adolescent is different and unknown, making cessation the only safe goal."
        },
        'IV': {
            "text": "Vaping’s risks and benefits remain poorly understood in children, but it has been shown vaping has clear benefits over cigarettes in children.",
            "evaluation": "INCORRECT. The claim of 'clear benefits' for children is false and dangerous. This approach is not endorsed by public health or medical bodies."
        },
        'V': {
            "text": "Consider initiating bupropion and varenicline depending on her son’s needs.",
            "evaluation": "CORRECT. These prescription medications are valid considerations for a physician treating adolescent nicotine dependence, often as a second-line or for more severe cases."
        }
    }

    print("Step 1: Evaluate each individual statement.")
    print("-" * 50)
    for key, value in statements.items():
        print(f"Statement {key}: {value['text']}")
        print(f"Evaluation: {value['evaluation']}\n")

    print("Step 2: Identify the best combination of counseling points.")
    print("-" * 50)
    print("Based on the evaluation, the most appropriate counseling points to consider are II, III, and V.")
    print("Options I and IV are inappropriate and potentially harmful.\n")
    print("We must now find the best combination from the available answer choices:")
    print("Choice J combines II and III.")
    print("  - III explains WHY the son must quit (risks to adolescents are different and unknown).")
    print("  - II provides a WHAT for the next step (use NRT as a first-line cessation tool).\n")
    print("This combination is the most comprehensive and effective for initial counseling. It provides both the crucial rationale and a concrete, actionable first step.\n")
    
    final_answer = 'J'
    print(f"Final Conclusion: The best answer choice is the one that combines statements II and III.")
    
    # Per instructions, the final output must be in the format <<<answer content>>>
    print("The final selected option corresponds to choice J.")

    # The actual answer to be parsed
    return f"<<<{final_answer}>>>"

# Execute the analysis and print the final answer
final_answer_formatted = analyze_vaping_counseling_options()
# The final line below is commented out to prevent double printing in the final output,
# but the function itself returns the required format. The 'print' calls inside the function handle the explanation.
# print(final_answer_formatted)
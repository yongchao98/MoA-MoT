import textwrap

def solve_case():
    """
    Analyzes the clinical case to determine the best next diagnostic step.
    """
    case_summary = """
    The patient presents with symptoms and signs highly characteristic of allergic contact dermatitis from clothing (textile dermatitis). The key features are the eczematous rash distribution on the periphery of the axillary folds (sparing the vaults) and the history of wearing new workout clothes.
    """

    question = "Which of the following is the best next step in diagnosis?"

    options = {
        'A': "Skin biopsy",
        'B': "KOH preparation",
        'C': "Topical steroid",
        'D': "Patch test",
        'E': "None of the above"
    }

    # Analysis of options
    analysis = {
        'A': "A skin biopsy is generally not the first-line diagnostic test for suspected contact dermatitis; it's more invasive.",
        'B': "A KOH preparation tests for fungal infections, which is less likely given the specific clinical picture.",
        'C': "A topical steroid is a form of treatment, not a step in diagnosis.",
        'D': "A patch test is the definitive, gold-standard method to identify the specific allergen causing allergic contact dermatitis. The case text itself confirms this was the procedure used to find the cause ('Patch testing was performed, and positive reactions were observed...')."
    }

    correct_option_key = 'D'
    correct_option_value = options[correct_option_key]
    reasoning = analysis[correct_option_key]

    print("Case Analysis:")
    print(textwrap.fill(case_summary.strip(), width=80))
    print("\n" + question)
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\nConclusion:")
    print(textwrap.fill(reasoning, width=80))
    print(f"\nTherefore, the best next step is D: {correct_option_value}.")

solve_case()

print("<<<D>>>")
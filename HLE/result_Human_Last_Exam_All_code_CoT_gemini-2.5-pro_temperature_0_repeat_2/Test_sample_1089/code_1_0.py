def solve_turbine_blade_question():
    """
    Analyzes the sources of damage to turbine blades and identifies which is
    best addressed by TIG welding build-up repair.
    """
    question = "What is the main source of damage addressed by manual TIG welding repair (build-up of layers of filler material)?"

    options = {
        'A': 'Stress Corrosion Cracking',
        'B': 'Foreign Object Damage',
        'C': 'Blade Tip Rub and Wear',
        'D': 'Creep Deformation',
        'E': 'Fatigue Cracking',
        'F': 'High-Temperature Oxidation and Corrosion'
    }

    analysis = {
        'A': "Incorrect. Cracking is a fracture, not material loss. While welding can be part of a crack repair (after grinding out the crack), the primary damage isn't addressed by 'building up' material.",
        'B': "Plausible. FOD creates nicks and gouges, which are forms of material loss. Welding is used to fill these. This is a valid application.",
        'C': "Most Correct. This is a classic and frequent repair. Blade tips lose material due to rubbing against the engine casing. TIG welding is the standard procedure to 'build-up' the tip, restoring its length and the required engine clearance. This perfectly matches the description.",
        'D': "Incorrect. Creep is the slow stretching and distortion of the entire blade under heat and stress. It is not a localized material loss and is generally not repairable by welding; the blade is replaced.",
        'E': "Incorrect. Similar to other forms of cracking, the primary damage is a fracture. The repair focuses on stopping the crack, not just building up lost material.",
        'F': "Plausible but less specific. This can cause surface pitting and thinning (material loss). However, the most common and direct application of 'building up layers' is for tip wear."
    }

    print("--- Analysis of Turbine Blade Damage Repair ---")
    print(f"Question: {question}\n")
    print("The key repair method is 'build-up of layers of filler material', which implies replacing lost material.\n")

    for key in options:
        print(f"Evaluating Option {key}: {options[key]}")
        print(f"   Analysis: {analysis[key]}\n")

    correct_answer_key = 'C'
    print("--- Final Conclusion ---")
    print(f"The source of damage that most directly and commonly corresponds to a 'build-up' repair is '{options[correct_answer_key]}'.")
    print("This process restores material lost at the blade tip, which is a routine maintenance procedure.")
    print(f"Therefore, the final answer is {correct_answer_key}.")


solve_turbine_blade_question()
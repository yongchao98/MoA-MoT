def solve_turbine_blade_question():
    """
    Analyzes the repair of aeroengine turbine blades to determine the primary
    damage source addressed by TIG welding build-up.
    """

    choices = {
        'A': "Stress Corrosion Cracking",
        'B': "Foreign Object Damage",
        'C': "Blade Tip Rub and Wear",
        'D': "Creep Deformation",
        'E': "Fatigue Cracking",
        'F': "High-Temperature Oxidation and Corrosion"
    }

    # The repair method described, "build-up of layers of filler material,"
    # is an additive process designed to replace lost material and restore geometry.
    
    # Analysis: While TIG welding can repair cracks (A, E) or damage from foreign objects (B),
    # a very common and routine application of material build-up is to address the wear
    # on the tip of the turbine blade. This wear (C) occurs from contact with the engine
    # casing and reduces engine efficiency. Restoring the blade tip length is a primary MRO task.
    
    correct_choice_key = 'C'
    explanation = (
        "The key is the method: 'build-up of layers'. This directly addresses material loss. "
        "Blade Tip Rub and Wear is a common, predictable operational issue that results in "
        "material loss at the blade tip. Restoring this lost material with TIG welding is a "
        "standard MRO procedure to maintain engine efficiency."
    )

    print("Explanation:", explanation)
    print("-" * 20)
    print("The final answer is:")
    
    # As requested, printing the components of the final answer.
    print(correct_choice_key)
    print(choices[correct_choice_key])

solve_turbine_blade_question()
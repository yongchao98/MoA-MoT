def solve_turbine_blade_problem():
    """
    Analyzes the relationship between TIG welding repair and types of turbine blade damage
    to determine the most likely application.
    """
    question = "What is the main source of damage addressed by manual TIG welding repair by build-up of layers?"

    damage_options = {
        'A': "Stress Corrosion Cracking",
        'B': "Foreign Object Damage",
        'C': "Blade Tip Rub and Wear",
        'D': "Creep Deformation",
        'E': "Fatigue Cracking",
        'F': "High-Temperature Oxidation and Corrosion"
    }

    # The repair method is "build-up of layers of filler material".
    # This implies restoring lost material to regain original dimensions.

    # Analysis:
    # - Cracking (A, E) is about filling a fracture, not large-scale material build-up.
    # - Creep (D) is a distortion/stretching of the entire blade, not typically fixed by adding material.
    # - Oxidation/Corrosion (F) is a surface issue, often treated with new coatings.
    # - Blade Tip Rub and Wear (C) is the literal grinding away of material from the blade tip.
    #   This lost material must be replaced to restore the blade's height and profile.
    #   This perfectly matches the "build-up of layers" repair description.

    correct_answer_key = 'C'
    explanation = (
        "The repair method described, 'build-up of layers of filler material' via TIG welding, "
        "is an additive process used to restore material that has been lost.\n\n"
        "Of the choices provided, 'Blade Tip Rub and Wear' is the most fitting application. This damage occurs "
        "when the rotating blade tips make contact with the stationary casing, progressively grinding away "
        "the tip material. To repair this, welders meticulously add, or 'build-up', layers of a "
        "compatible alloy onto the blade tip to restore its original height. The blade is then machined "
        "back to its precise aerodynamic profile. This common MRO procedure directly addresses material loss and "
        "restores the engine's performance and efficiency."
    )

    print(f"Question: {question}")
    print("-" * 30)
    print(f"Correct Answer Choice: {correct_answer_key}")
    print(f"Which corresponds to: {damage_options[correct_answer_key]}")
    print("-" * 30)
    print("Detailed Explanation:")
    print(explanation)

solve_turbine_blade_problem()
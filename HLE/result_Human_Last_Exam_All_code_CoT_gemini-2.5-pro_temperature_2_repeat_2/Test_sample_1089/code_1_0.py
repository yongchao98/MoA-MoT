import textwrap

def solve_turbine_blade_question():
    """
    Analyzes the repair process of aeroengine turbine blades to determine the primary
    damage source addressed by TIG welding build-up.
    """

    question = "What is the main source of damage addressed by manual TIG welding repair involving the build-up of layers of filler material on aeroengine turbine blades?"

    options = {
        'A': "Stress Corrosion Cracking",
        'B': "Foreign Object Damage",
        'C': "Blade Tip Rub and Wear",
        'D': "Creep Deformation",
        'E': "Fatigue Cracking",
        'F': "High-Temperature Oxidation and Corrosion"
    }

    # Analysis of the repair method
    # The key term is "build-up of layers." This implies a process of restoring volume and
    # geometry to a component where material has been physically removed or worn away,
    # rather than just joining a crack.

    # Evaluating the options based on the analysis
    # - Cracking (A, E) is repaired by welding, but the primary action is fusing a line, not a bulk build-up.
    # - Creep (D) is a bulk deformation and is not typically fixed by welding.
    # - FOD (B) and Corrosion (F) can cause material loss requiring build-up.
    # - However, Blade Tip Rub and Wear (C) is a highly common and predictable type of damage
    #   where material is systematically ground off the blade tip. The standard repair is to
    #   weld new material onto the tip to restore its original length and the critical
    #   tip-to-casing clearance. This procedure is a quintessential example of "build-up."

    correct_answer_key = 'C'
    correct_answer_text = options[correct_answer_key]

    explanation = (
        "The repair process described as 'build-up of layers of filler material' "
        "is primarily used to restore material that has been lost from a surface. "
        "Blade Tip Rub and Wear is a common mechanism in aeroengines where the tip of the "
        "turbine blade makes contact with the casing, grinding away material. "
        "The standard repair procedure for this is to use manual TIG welding to "
        "add material back onto the tip, thereby 'building it up' to its original dimensions."
    )

    print("Question:", question)
    print("-" * 30)
    print("Analysis:")
    print("\n".join(textwrap.wrap(explanation, width=70)))
    print("-" * 30)
    print(f"The most fitting answer is:\n({correct_answer_key}) {correct_answer_text}")

if __name__ == '__main__':
    solve_turbine_blade_question()
import textwrap

def solve_turbine_blade_question():
    """
    Analyzes the repair method for aeroengine turbine blades to determine the
    primary source of damage it addresses.
    """
    question = "What is the main source of damage addressed by manual TIG welding (GTAW) build-up of layers of filler material?"

    choices = {
        'A': "Stress Corrosion Cracking",
        'B': "Foreign Object Damage",
        'C': "Blade Tip Rub and Wear",
        'D': "Creep Deformation",
        'E': "Fatigue Cracking",
        'F': "High-Temperature Oxidation and Corrosion"
    }

    # Analysis: The key is "build-up of layers of filler material".
    # This implies restoring geometry or replacing material that has been worn away.

    reasoning = {
        'A': "Incorrect. Repairing cracks typically involves removing the crack and filling the void, not a general build-up of layers.",
        'B': "Plausible, but not primary. FOD repair involves filling localized nicks and gouges.",
        'C': "Correct. Blade tip wear is a common issue where material is lost from the tip. The standard repair is to add layers of weld material to build the tip back to its original dimensions, restoring critical engine clearances.",
        'D': "Incorrect. Creep is a bulk deformation of the material that welding cannot reverse. Crept blades are often retired.",
        'E': "Incorrect. Same reasoning as for Stress Corrosion Cracking (A).",
        'F': "Plausible, but not primary. While build-up can address material loss from corrosion, it's the defining repair for predictable and performance-critical tip wear."
    }

    # The most direct and common application of "building up layers" is to restore blade tips.
    correct_answer_key = 'C'

    print("Step-by-step analysis of the problem:")
    print("-" * 70)
    print("The specified repair method is 'manual TIG welding build-up'. This technique is used to add material to a component.")
    print("We need to find the damage type for which this is the primary repair method.")
    print("\nEvaluation of options:")
    for key, text in choices.items():
        print(f"  - Option {key} ({text}):")
        # Wrap the reasoning text for better readability
        wrapped_text = textwrap.fill(reasoning[key], width=65, initial_indent='    ', subsequent_indent='    ')
        print(wrapped_text)

    print("-" * 70)
    print("Conclusion: The most characteristic and primary use of TIG build-up on turbine blades is to restore material lost due to tip rub and wear.")
    print("\nFinal Answer Equation:")
    print(f"{correct_answer_key} = {choices[correct_answer_key]}")


solve_turbine_blade_question()
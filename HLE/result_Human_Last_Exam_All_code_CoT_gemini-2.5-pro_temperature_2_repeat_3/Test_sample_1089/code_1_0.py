import textwrap

def solve_turbine_blade_problem():
    """
    Analyzes the provided question about turbine blade repair and determines the most likely answer.
    """
    question = "What is the main source of damage addressed by manual TIG welding repair involving build-up of layers of filler material?"

    options = {
        'A': "Stress Corrosion Cracking",
        'B': "Foreign Object Damage",
        'C': "Blade Tip Rub and Wear",
        'D': "Creep Deformation",
        'E': "Fatigue Cracking",
        'F': "High-Temperature Oxidation and Corrosion"
    }

    # Analysis of the repair method
    analysis_repair_method = textwrap.fill(
        "The repair method described is 'manual TIG welding (GTAW) build-up of layers of filler material'. "
        "This explicitly points to a process where lost material is added back to restore a component's original dimensions.",
        width=100
    )

    # Evaluation of each damage type
    evaluations = {
        'A': "Repairing cracks might use welding, but 'build-up of layers' better describes restoring worn volume.",
        'B': "Possible, but tip rub is a more regular and predictable wear mechanism requiring dimensional restoration in turbines.",
        'C': "This is the classic case. Blade tips wear down from contact with the engine casing, losing material. TIG welding is the standard MRO procedure to 'build-up' the tip to its original height, restoring engine efficiency. This is a perfect match.",
        'D': "Creep is plastic stretching of the entire blade; it cannot be fixed by adding material to one spot with welding.",
        'E': "Similar to stress cracking, the description points more to restoring lost volume than fixing a linear crack.",
        'F': "Usually addressed by reapplying protective coatings. Only severe localized pitting would require a weld build-up."
    }

    # Conclusion
    correct_answer_key = 'C'
    conclusion = textwrap.fill(
        f"Based on the analysis, {options[correct_answer_key]} is the most direct and common reason for performing a TIG welding 'build-up' repair on turbine blades. "
        "The process perfectly addresses the loss of material from the blade tip.",
        width=100
    )

    print("Problem Analysis:")
    print("-" * 20)
    print(analysis_repair_method)
    print("\nEvaluation of Answer Choices:")
    print("-" * 20)
    for key, text in evaluations.items():
        print(f"({key}) {options[key]}: {text}")

    print("\nFinal Conclusion:")
    print("-" * 20)
    print(conclusion)
    print(f"\nTherefore, the correct answer choice is {correct_answer_key}.")


solve_turbine_blade_problem()
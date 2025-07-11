def solve_agitation_case():
    """
    Analyzes treatment options for an agitated patient using a scoring model.
    The scenario: A violent patient in the ER, with an unknown history,
    is unresponsive to an initial 5mg IM olanzapine (Zyprexa) dose.
    """

    options = {
        "A": {
            "description": "2mg IV lorazepam",
            "scores": {"Efficacy": 1, "Patient Safety": -1, "Staff Safety": -2, "Appropriateness": 0},
            "reasoning": "IV access is dangerous and difficult in a violent patient, posing a high risk of needlestick injury to staff. IV benzodiazepines also carry a higher risk of respiratory depression."
        },
        "B": {
            "description": "2mg IM lorazepam + 5mg olanzapine IM",
            "scores": {"Efficacy": 2, "Patient Safety": 1, "Staff Safety": 2, "Appropriateness": 2},
            "reasoning": "This is an excellent choice. It escalates care appropriately by adding a synergistic agent (lorazepam) and increasing the olanzapine to a full 10mg total dose. The IM route is safe for staff, and this combination is highly effective for refractory agitation."
        },
        "C": {
            "description": "Verbal de-escalation before any pharmacologic intervention",
            "scores": {"Efficacy": -2, "Patient Safety": 0, "Staff Safety": -1, "Appropriateness": -2},
            "reasoning": "The patient has already escalated to physical violence. The window for verbal de-escalation as a primary strategy has passed. Prioritizing it now would endanger staff."
        },
        "D": {
            "description": "10mg IM olanzapine",
            "scores": {"Efficacy": 0, "Patient Safety": 2, "Staff Safety": 2, "Appropriateness": 1},
            "reasoning": "This involves increasing the dose of a medication that has already proven ineffective at a lower dose. While safe, it is less likely to be effective than combination therapy."
        },
        "E": {
            "description": "10mg IM olanzapine + 2mg IM lorazepam",
            "scores": {"Efficacy": 2, "Patient Safety": -1, "Staff Safety": 2, "Appropriateness": 0},
            "reasoning": "This option implies giving an *additional* 10mg of olanzapine (total 15mg), which is an aggressive and potentially unsafe dose in a patient with an unknown medical history."
        }
    }

    best_option = None
    max_score = -float('inf')

    print("Evaluating options based on a clinical scoring model (-2=Very Poor, +2=Excellent)...\n")

    for key, data in options.items():
        scores = data["scores"]
        total_score = sum(scores.values())

        # Print the "equation" for each option
        print(f"Option {key}: {data['description']}")
        print(f"   Scoring Equation: Efficacy({scores['Efficacy']}) + Patient Safety({scores['Patient Safety']}) + Staff Safety({scores['Staff Safety']}) + Appropriateness({scores['Appropriateness']}) = {total_score}")
        print("-" * 20)

        if total_score > max_score:
            max_score = total_score
            best_option = key

    print(f"\nConclusion:")
    print(f"The best next step is Option {best_option} with a score of {max_score}.")
    print(f"Reasoning: {options[best_option]['reasoning']}")


solve_agitation_case()
<<<B>>>
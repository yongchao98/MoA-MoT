import operator

def evaluate_agitation_options():
    """
    Evaluates and scores different treatment options for an agitated patient
    who failed initial monotherapy, based on clinical criteria.
    """
    options = {
        "A. 2mg IV lorazepam": {
            "Efficacy": 4,
            "Safety": 2,
            "Feasibility": 1
        },
        "B. 2mg IM lorazepam + 5mg olanzapine IM": {
            "Efficacy": 5,
            "Safety": 3,
            "Feasibility": 4
        },
        "C. Verbal de-escalation before any pharmacologic intervention": {
            "Efficacy": 1,
            "Safety": 5,
            "Feasibility": 1
        },
        "D. 10mg IM olanzapine": {
            "Efficacy": 3,
            "Safety": 4,
            "Feasibility": 4
        },
        "E. 10mg IM olanzapine + 2mg IM lorazepam": {
            "Efficacy": 5,
            "Safety": 1,
            "Feasibility": 4
        }
    }

    scores = {}
    print("Scoring each option based on Efficacy, Safety, and Feasibility (scale 1-5):\n")

    for option, criteria in options.items():
        efficacy = criteria["Efficacy"]
        safety = criteria["Safety"]
        feasibility = criteria["Feasibility"]
        total_score = efficacy + safety + feasibility
        scores[option] = total_score

        # Output each number in the final equation
        print(f"Scoring for {option}")
        print(f"Equation: Efficacy({efficacy}) + Safety({safety}) + Feasibility({feasibility}) = {total_score}\n")

    # Determine the best option
    best_option = max(scores.items(), key=operator.itemgetter(1))

    print("--------------------------------------------------")
    print(f"The best next step is the option with the highest score: '{best_option[0]}'")
    print(f"This approach combines high efficacy with acceptable safety and feasibility after initial treatment has failed.")
    print("--------------------------------------------------")


if __name__ == "__main__":
    evaluate_agitation_options()
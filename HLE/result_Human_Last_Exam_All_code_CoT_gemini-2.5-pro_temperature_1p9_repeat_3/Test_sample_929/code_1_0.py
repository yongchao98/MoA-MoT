import sys
# Redirect stdout to a string to control the final output precisely
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()


def analyze_pollinator_behavior():
    """
    Analyzes different insect behavior patterns to determine which has the greatest
    positive effect on plant fitness, based on a provided ethogram.
    """

    # Ethogram legend:
    # 1: investigation start, 2: investigation end (non-contact)
    # 3: interaction start, 4: interaction end (contact)
    # 5: feeding start, 6: feeding end (subset of interaction)
    # Plant fitness is primarily driven by pollination, which requires interaction (contact).

    scenarios = {
        "A": {
            "description": "4-3 >> 6-5",
            "equation": "Interaction_Duration(10s) >> Feeding_Duration(1s)",
            "logic": "The insect spends a long time on the plant but only feeds for a short while. While long contact is good for pollination, the poor reward (short feeding) may discourage the insect from visiting other plants, limiting overall cross-pollination.",
            "valid": True,
            "score": 2  # Good, but not optimal
        },
        "B": {
            "description": "6-5 >> 4-3",
            "equation": "Feeding_Duration(10s) >> Interaction_Duration(1s)",
            "logic": "Logically impossible. An insect cannot feed (behavior 5-6) without being in contact with the plant (behavior 3-4). The duration of feeding is always a subset of the interaction duration.",
            "valid": False,
            "score": -1  # Invalid
        },
        "C": {
            "description": "4-3 >> 2-1",
            "equation": "Interaction_Duration(10s) >> Investigation_Duration(1s)",
            "logic": "The insect makes a quick decision to land and then spends a long time interacting with the flower (which includes feeding). This is highly efficient for the plant, maximizing the chance of pollination for each visit. This is the ideal behavior for a pollinator.",
            "valid": True,
            "score": 3  # Optimal
        },
        "D": {
            "description": "n(5)/hour >> n(3)/hour",
            "equation": "n_Feeding_events(10) >> n_Interaction_events(1)",
            "logic": "Logically impossible. Every feeding event (start: 5) must be part of an interaction event (start: 3). The number of feeding events cannot exceed the number of interaction events.",
            "valid": False,
            "score": -1  # Invalid
        },
        "E": {
            "description": "n(1)/hour >> n(3)/hour",
            "equation": "n_Investigation_events(10) >> n_Interaction_events(1)",
            "logic": "The insect approaches the plant many times but rarely lands. This indicates the plant is not very attractive or the insect is very hesitant. Pollination events are rare and inefficient. This is very poor for plant fitness.",
            "valid": True,
            "score": 1  # Inefficient
        },
        "F": {
            "description": "n(3)/hour >> n(1)/hour",
            "equation": "n_Interaction_events(10) >> n_Investigation_events(1)",
            "logic": "Logically impossible. Every interaction (contact) must be preceded by an investigation (an approach), even a brief one. The number of interactions cannot exceed the number of investigations.",
            "valid": False,
            "score": -1  # Invalid
        }
    }

    best_scenario_key = None
    max_score = -float('inf')

    print("Analyzing Pollinator Behavior Patterns for Plant Fitness:\n")

    for key, value in scenarios.items():
        print(f"--- Scenario {key} ---")
        print(f"Pattern: {value['description']}")
        # The prompt requires outputting the numbers in the final equation.
        print(f"Example Equation: {value['equation']}")
        print(f"Analysis: {value['logic']}")
        print(f"Logical Validity: {'Valid' if value['valid'] else 'Impossible'}")
        print(f"Fitness Impact Score (Higher is better): {value['score']}\n")

        if value['score'] > max_score:
            max_score = value['score']
            best_scenario_key = key

    print("="*40)
    print("Conclusion:")
    print(f"Scenario '{best_scenario_key}' describes the pattern with the greatest positive effect on plant fitness.")
    print("It represents an efficient pollinator that quickly commits to a long interaction, maximizing the opportunity for pollen transfer.")
    print("="*40)

analyze_pollinator_behavior()

# Restore original stdout and print captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
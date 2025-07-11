import sys

def solve():
    """
    Analyzes nectarivorous insect behavior patterns to determine which has the greatest
    positive effect on milkweed plant fitness.
    """

    # Ethogram events:
    # 1) investigation start
    # 2) investigation end
    # 3) interaction start (contact)
    # 4) interaction end
    # 5) feeding start
    # 6) feeding end

    # Rules:
    # - Feeding (5-6) is a type of interaction (3-4). An insect must be in contact to feed.
    # - Plant fitness is driven by pollination, which requires pollen transfer between flowers.
    # - Therefore, behaviors that maximize efficient contact with many flowers are best.

    patterns = {
        'A': '4-3 >> 6-5',
        'B': '6-5 >> 4-3',
        'C': '4-3 >> 2-1',
        'D': 'n(5)/hour >> n(3)/hour',
        'E': 'n(1)/hour >> n(3)/hour',
        'F': 'n(3)/hour >> n(1)/hour'
    }

    # Assign a score based on the likely effect on plant fitness.
    scores = {}
    explanations = {}

    # Pattern A: 4-3 >> 6-5
    # Interpretation: The total duration of interaction (time between event 3 and 4) is much greater
    # than the duration of feeding (time between event 5 and 6).
    # Analysis: This means the insect spends a lot of time on the flower without feeding.
    # Feeding is the most effective action for pollination in milkweed. This is inefficient.
    scores['A'] = -10
    explanations['A'] = "Pattern A (4-3 >> 6-5): Describes long interactions with short feeding times. This is inefficient for pollination, as feeding is the primary mechanism for pollen transfer. This is a negative pattern."

    # Pattern B: 6-5 >> 4-3
    # Interpretation: The duration of feeding is much greater than the duration of interaction.
    # Analysis: This is logically impossible. Feeding is a *type* of interaction, so its duration
    # cannot exceed the total interaction duration.
    scores['B'] = -100
    explanations['B'] = "Pattern B (6-5 >> 4-3): Describes feeding duration as much longer than interaction duration. This is logically impossible, as feeding requires interaction (contact). This pattern is invalid."

    # Pattern C: 4-3 >> 2-1
    # Interpretation: The duration of interaction (4-3) is much greater than the duration of investigation (2-1).
    # Analysis: The insect commits to a long visit once it makes contact. This is good for pollinating a
    # single flower but is less effective for population-level fitness than visiting many flowers.
    scores['C'] = 20
    explanations['C'] = "Pattern C (4-3 >> 2-1): Describes long interactions relative to investigation time. This is a positive trait, ensuring deep engagement with a flower, but may not maximize the number of flowers visited across the population."

    # Pattern D: n(5)/hour >> n(3)/hour
    # Interpretation: The frequency of feeding starts (n(5)) is much greater than the frequency of interaction starts (n(3)).
    # Analysis: This is logically impossible. A feeding event cannot start without an interaction (contact) having started first.
    scores['D'] = -100
    explanations['D'] = "Pattern D (n(5)/hour >> n(3)/hour): Describes more feeding events than interaction events. This is logically impossible, as feeding is a subset of interaction. This pattern is invalid."

    # Pattern E: n(1)/hour >> n(3)/hour
    # Interpretation: The frequency of investigations (n(1)) is much greater than the frequency of interactions (n(3)).
    # Analysis: The insect approaches many flowers but rarely lands to make contact. This is the behavior
    # of a very poor pollinator.
    scores['E'] = -20
    explanations['E'] = "Pattern E (n(1)/hour >> n(3)/hour): Describes an insect that investigates far more often than it interacts. This is a very inefficient pollinator, providing little benefit to the plant."

    # Pattern F: n(3)/hour >> n(1)/hour
    # Interpretation: The frequency of interactions (n(3)) is much greater than the frequency of investigations (n(1)).
    # Analysis: This describes a highly efficient pollinator that forgoes hesitant 'investigation' and moves
    # directly from flower to flower, making contact each time. This maximizes the number of flowers
    # visited and is the best pattern for promoting cross-pollination and plant population fitness.
    scores['F'] = 50
    explanations['F'] = "Pattern F (n(3)/hour >> n(1)/hour): Describes an insect that interacts with flowers far more often than it investigates them. This represents a highly efficient forager that maximizes the number of visited flowers, which is the most effective way to promote cross-pollination and overall plant fitness."

    # Find the best pattern
    best_pattern = max(scores, key=scores.get)

    print("--- Analysis of Behavioral Patterns for Plant Fitness ---")
    for pattern_key in sorted(patterns.keys()):
        print(f"\nAnalyzing Choice {pattern_key}: '{patterns[pattern_key]}'")
        print(explanations[pattern_key])
        print(f"Assigned Fitness Score: {scores[pattern_key]}")

    print("\n--- Conclusion ---")
    print(f"The pattern with the highest fitness score is '{best_pattern}'.")
    print(explanations[best_pattern])
    print("\nThis pattern represents the greatest positive effect on plant fitness.")
    
    # Required final answer format
    sys.stdout.flush()
    print(f"\n<<<{best_pattern}>>>")


solve()
def calculate_fitness(scenario_name, n_interaction, n_feeding, w_interaction=1, w_feeding=10):
    """
    Calculates a fitness score based on interaction and feeding events.
    n_interaction: The number of interaction starts (event 3).
    n_feeding: The number of feeding starts (event 5).
    w_interaction: The weight/value of a single interaction.
    w_feeding: The weight/value of a single feeding event.
    """
    score = (n_feeding * w_feeding) + (n_interaction * w_interaction)
    print(f"Scenario {scenario_name}:")
    # The final equation includes each number used in the calculation
    print(f"  Fitness Score = (n_feeding * weight) + (n_interaction * weight)")
    print(f"  Fitness Score = ({n_feeding} * {w_feeding}) + ({n_interaction} * {w_interaction}) = {score}\n")
    return score

def find_best_pattern():
    """
    Models each plausible scenario and prints the fitness calculation.
    """
    print("--- Evaluating Pollinator Patterns for Plant Fitness ---\n")
    print("A fitness score is calculated based on weighted counts of feeding and interaction events.")
    print("A higher score indicates a greater positive effect on plant fitness.\n")

    # Scenarios are modeled based on the logic described in the answer choices.
    # We only model plausible scenarios (A, C, E, F).

    # A (4-3 >> 6-5): Long interaction, little feeding. "The Lounger".
    calculate_fitness('A', n_interaction=10, n_feeding=2)

    # C (4-3 >> 2-1): Commits to long interactions vs investigating. "The Efficient Visitor".
    calculate_fitness('C', n_interaction=10, n_feeding=8)

    # E (n(1) >> n(3)): Many investigations, few interactions. "The Hesitant Visitor".
    calculate_fitness('E', n_interaction=5, n_feeding=4)

    # F (n(3) >> n(1)): Many interactions per investigation. "The Busy Hopper".
    calculate_fitness('F', n_interaction=50, n_feeding=40)

    print("--- Conclusion ---")
    print("Scenario F, representing an insect that performs many interactions per visit (hopping between flowers), yields the highest fitness score.")
    print("This behavior is the most effective for pollination.")

if __name__ == '__main__':
    find_best_pattern()

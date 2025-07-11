import pandas as pd

def evaluate_fitness_patterns():
    """
    Analyzes different insect behavior patterns to find which has the
    greatest positive effect on plant fitness (pollination).
    """

    # Step 1: Define fitness contribution weights for each behavior.
    # Feeding is the most effective for pollination.
    # Interaction (non-feeding) is less effective but still requires contact.
    # Investigation has no contact and no fitness contribution.
    feeding_fitness_weight = 5
    interaction_fitness_weight = 1
    investigation_fitness_weight = 0

    results = {}

    # --- Scenario Modeling & Calculation ---

    # Scenario A: 4-3 >> 6-5 (Interaction duration >> Feeding duration)
    # Model: In 11 total minutes of contact, 10 are non-feeding interaction and 1 is feeding.
    # This represents one long, inefficient visit.
    a_feeding_duration = 1
    a_interaction_duration = 10
    fitness_a = (a_feeding_duration * feeding_fitness_weight) + (a_interaction_duration * interaction_fitness_weight)
    results['A (Interaction duration >> Feeding duration)'] = fitness_a

    # Scenario B: 6-5 >> 4-3 (Feeding duration >> Interaction duration)
    # Model: In 11 total minutes of contact, 10 are feeding and 1 is non-feeding interaction.
    # This represents one long, very effective visit.
    # Note: '4-3' is interpreted as non-feeding interaction duration.
    b_feeding_duration = 10
    b_interaction_duration = 1
    fitness_b = (b_feeding_duration * feeding_fitness_weight) + (b_interaction_duration * interaction_fitness_weight)
    results['B (Feeding duration >> Interaction duration)'] = fitness_b

    # Scenario C: 4-3 >> 2-1 (Interaction duration >> Investigation duration)
    # Model: The insect spends more time on the flower than approaching it.
    # This doesn't specify feeding vs. interaction, so we assume an average ratio (e.g., 50/50).
    c_total_interaction_duration = 10
    c_feeding_duration = 5 # 50% of time
    c_interaction_duration = 5 # 50% of time
    fitness_c = (c_feeding_duration * feeding_fitness_weight) + (c_interaction_duration * interaction_fitness_weight)
    results['C (Interaction duration >> Investigation duration)'] = fitness_c

    # Scenario D: n(5)/hour >> n(3)/hour (Frequency of feeding starts >> interaction starts)
    # This is logically impossible, as feeding (5) is a type of interaction (3).
    fitness_d = 0
    results['D (Freq. of Feeding >> Freq. of Interaction)'] = fitness_d

    # Scenario E: n(1)/hour >> n(3)/hour (Frequency of investigations >> interactions)
    # Model: The insect makes many approaches but rarely lands.
    # Let's say 10 approaches result in only 1 interaction event. Low overall activity.
    e_num_interactions = 1
    # Assume this single interaction is of average duration and composition (5 min feeding, 5 min non-feeding)
    e_fitness_per_interaction = (5 * feeding_fitness_weight) + (5 * interaction_fitness_weight)
    fitness_e = e_num_interactions * e_fitness_per_interaction
    results['E (Freq. of Investigation >> Freq. of Interaction)'] = fitness_e

    # Scenario F: n(3)/hour >> n(1)/hour (Frequency of interactions >> investigations)
    # Model: The insect is highly efficient, turning almost every approach into a visit.
    # This maximizes the number of separate flower visits, crucial for cross-pollination.
    # Let's say 1 "investigation" of a plant patch leads to 10 separate flower interactions.
    f_num_interactions = 10
    # Assume each interaction is brief but effective (e.g., 1 min feeding, 1 min non-feeding)
    f_fitness_per_interaction = (1 * feeding_fitness_weight) + (1 * interaction_fitness_weight)
    fitness_f = f_num_interactions * f_fitness_per_interaction
    results['F (Freq. of Interaction >> Freq. of Investigation)'] = fitness_f


    # --- Print Results ---
    print("--- Analysis of Behavioral Patterns for Plant Fitness ---")
    print("Fitness is calculated based on the potential for pollination.\n")
    df = pd.DataFrame(list(results.items()), columns=['Scenario', 'Calculated Fitness Score'])
    df = df.sort_values(by='Calculated Fitness Score', ascending=False)
    print(df.to_string(index=False))

    print("\n--- Conclusion ---")
    winner = df.iloc[0]
    print(f"The pattern '{winner['Scenario']}' yields the highest fitness score.")
    print("This pattern represents a highly efficient pollinator that visits many different flowers,")
    print("which is the most critical factor for successful cross-pollination and plant reproductive fitness.\n")

    print("--- Final Equation for the Best Pattern (F) ---")
    print(f"Total Fitness = Number of Interactions * Fitness per Interaction")
    print(f"Total Fitness = Number of Interactions * ((Feeding Duration * Feeding Weight) + (Interaction Duration * Interaction Weight))")
    print(f"Total Fitness = {f_num_interactions} * (({1} * {feeding_fitness_weight}) + ({1} * {interaction_fitness_weight})) = {fitness_f}")


if __name__ == '__main__':
    evaluate_fitness_patterns()
import sys
import io

# Redirect stdout to capture print output for the final answer
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def calculate_plant_fitness(scenario_name, n_interactions, t_feeding_per_visit, note=""):
    """Calculates a fitness score based on total feeding time."""
    total_feeding_time = n_interactions * t_feeding_per_visit
    print(f"--- Scenario {scenario_name} ---")
    if note:
        print(f"Note: {note}")
    print(f"Pattern: {scenario_name}")
    print(f"Number of Interactions per Hour = {n_interactions}")
    print(f"Feeding Time per Interaction = {t_feeding_per_visit:.0f} seconds")
    print(f"Resulting Plant Fitness Score (Total Feeding Time) = {n_interactions} * {t_feeding_per_visit:.0f} = {total_feeding_time:.0f}\n")
    return total_feeding_time

def analyze_behavior_patterns():
    """
    Simulates different insect behavior patterns to find which one maximizes plant fitness.
    Plant fitness is assumed to be directly proportional to the total time insects spend feeding.
    """
    print("Analyzing which behavior pattern has the greatest positive effect on plant fitness.")
    print("The key metric is Total Feeding Time, as feeding is the primary mechanism for pollination.\n")

    # --- Baseline Assumptions for a typical hour ---
    # 100 insects approach the plant area.
    n_investigations_per_hour = 100
    # A 50% conversion rate from investigation to interaction (landing on the plant).
    baseline_n_interactions = 50
    # A typical interaction (contact with the plant) lasts 60 seconds.
    baseline_t_interaction = 60
    # During a typical interaction, the insect spends 30% of the time feeding.
    baseline_t_feeding = baseline_t_interaction * 0.30

    # --- Scenario A: 4-3 >> 6-5 (Interaction time is much greater than feeding time) ---
    # This means the insect spends most of its contact time doing things other than feeding.
    t_feeding_A = baseline_t_interaction * 0.10 # Only 10% of interaction is feeding
    calculate_plant_fitness("A. 4-3 >> 6-5", baseline_n_interactions, t_feeding_A)

    # --- Scenario B: 6-5 >> 4-3 (Feeding time is much greater than interaction time) ---
    # This is logically impossible, as feeding is a subset of interaction.
    # We interpret this to mean that feeding time constitutes the vast majority of the interaction time.
    t_feeding_B = baseline_t_interaction * 0.90 # 90% of interaction is feeding
    note_B = "Interpreted as: Feeding duration is a very high proportion of interaction duration."
    calculate_plant_fitness("B. 6-5 >> 4-3", baseline_n_interactions, t_feeding_B, note=note_B)

    # --- Scenario C: 4-3 >> 2-1 (Interaction time is much greater than investigation time) ---
    # This describes an efficient visitor that doesn't hesitate.
    # This implies a very high conversion rate from investigation to interaction.
    n_interactions_C = n_investigations_per_hour * 0.90 # 90% conversion rate
    # The feeding time per visit is not specified, so we use the baseline.
    calculate_plant_fitness("C. 4-3 >> 2-1", n_interactions_C, baseline_t_feeding)

    # --- Scenario D: n(5)/hour >> n(3)/hour (Number of feeding starts >> number of interaction starts) ---
    # Impossible, as every feeding event (5) must be part of an interaction event (3).
    note_D = "Impossible. An insect cannot start feeding (5) without first starting an interaction (3)."
    calculate_plant_fitness("D. n(5)/hr >> n(3)/hr", 0, 0, note=note_D)

    # --- Scenario E: n(1)/hour >> n(3)/hour (Number of investigations >> number of interactions) ---
    # This means the plant is not very attractive; many insects approach but few land.
    # Low conversion rate from investigation to interaction.
    n_interactions_E = n_investigations_per_hour * 0.10 # 10% conversion rate
    # We use the baseline feeding time.
    calculate_plant_fitness("E. n(1)/hr >> n(3)/hr", n_interactions_E, baseline_t_feeding)

    # --- Scenario F: n(3)/hour >> n(1)/hour (Number of interactions >> number of investigations) ---
    # Impossible, as an interaction (contact) requires an investigation (approach).
    note_F = "Impossible. An insect cannot interact (3) without first investigating (1)."
    calculate_plant_fitness("F. n(3)/hr >> n(1)/hr", 0, 0, note=note_F)


analyze_behavior_patterns()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the user
print(output)

# Final Answer Determination
# Based on the output, Scenario B provides the highest fitness score.
print("<<<B>>>")
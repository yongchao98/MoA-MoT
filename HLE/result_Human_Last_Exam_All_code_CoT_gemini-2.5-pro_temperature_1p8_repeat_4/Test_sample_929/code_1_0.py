def calculate_and_print_fitness(scenario_name, n_interactions, duration_interaction, duration_feeding):
    """
    Calculates and prints a fitness score based on insect behavior, showing the equation.
    
    Assumption: Fitness is proportional to pollination. Feeding is 10x more effective
    for pollination than general, non-feeding contact.
    """
    # Pollination effectiveness constants (arbitrary units)
    k_feed = 10.0
    k_contact = 1.0

    # Ensure feeding duration isn't longer than interaction duration
    if duration_feeding > duration_interaction:
        print(f"Scenario '{scenario_name}' is logically impossible (feeding > interaction). Fitness = 0")
        return

    duration_non_feeding_contact = duration_interaction - duration_feeding
    
    # Pollination per single interaction event
    pollination_per_interaction = (k_contact * duration_non_feeding_contact) + (k_feed * duration_feeding)
    
    # Total pollination for all interactions in the period
    total_pollination = n_interactions * pollination_per_interaction
    
    print(f"--- {scenario_name} ---")
    print(f"Equation: Fitness = n_interactions * (k_contact * (d_interaction - d_feeding) + k_feed * d_feeding)")
    # Output each number in the final equation
    print(f"Fitness = {n_interactions} * ({k_contact} * ({duration_interaction} - {duration_feeding}) + {k_feed} * {duration_feeding}) = {total_pollination:.0f}")
    print("-" * 25)

# --- Main Program ---
print("Calculating plant fitness for different insect behavior patterns.")
print("The best pattern will have the highest fitness score.\n")

# Scenario C: 4-3 >> 2-1 (Long interaction duration >> short investigation duration)
# This implies an efficient visitor: high number of interactions and long duration.
# For example, 8 interactions of 120s each. Let's assume feeding is half the time.
calculate_and_print_fitness(
    scenario_name="C: Long Interaction >> Short Investigation",
    n_interactions=8,
    duration_interaction=120,
    duration_feeding=60
)

# Scenario A: 4-3 >> 6-5 (Long interaction duration >> short feeding duration)
# The insect stays on the plant but isn't performing the most useful action.
# Same parameters as C, but with very short feeding time.
calculate_and_print_fitness(
    scenario_name="A: Long Interaction >> Short Feeding",
    n_interactions=8,
    duration_interaction=120,
    duration_feeding=10
)

# Scenario E: n(1) >> n(3) (Many investigations, few interactions)
# An inefficient visitor that rarely makes contact.
# For example, only 2 interactions happen, with a standard duration.
calculate_and_print_fitness(
    scenario_name="E: Many Investigations >> Few Interactions",
    n_interactions=2,
    duration_interaction=60,
    duration_feeding=30
)

print("\nConclusion: Based on the fitness scores, Scenario C, where interaction time is")
print("much greater than investigation time, is the most beneficial for the plant.")
<<<C>>>
import math

def simulate_behavior_and_calculate_fitness():
    """
    Simulates different insect behavioral patterns over one hour
    to see which one has the greatest positive effect on plant fitness.
    """
    total_time_seconds = 3600  # 1 hour
    
    # Fitness model:
    # A successful interaction (contact) is good for pollination.
    # A feeding bout is even better, as it's the primary mechanism for pollen transfer.
    # Fitness Score = (num_interactions * 1) + (num_feedings * 5)
    
    print("Analyzing insect behavior patterns to determine the greatest positive effect on plant fitness.")
    print(f"Fitness will be scored over {total_time_seconds} seconds (1 hour).\n")

    scenarios = {
        'F': {
            'description': 'High efficiency: n(3)/hour >> n(1)/hour. Most investigations lead to interactions.',
            'cycle_time_s': 45,  # short visit: 5s invest, 10s interact, 30s travel
            'prob_interact': 0.9, # High conversion from investigation to interaction
            'prob_feed': 0.8,     # High likelihood of feeding during interaction
            'equation_template': "n(3)/hour ({:.1f}) >> n(1)/hour ({:.1f})"
        },
        'E': {
            'description': 'Low efficiency (Picky): n(1)/hour >> n(3)/hour. Few investigations lead to interactions.',
            'cycle_time_s': 45, # same short visit time
            'prob_interact': 0.1, # Low conversion from investigation to interaction
            'prob_feed': 0.8,
            'equation_template': "n(1)/hour ({:.1f}) >> n(3)/hour ({:.1f})"
        },
        'C': {
            'description': 'Long Interactions: 4-3 >> 2-1. Long time interacting vs. investigating.',
            'cycle_time_s': 155, # long visit: 5s invest, 120s interact, 30s travel
            'prob_interact': 0.9, # Not picky, just stays a long time
            'prob_feed': 0.8,
            'equation_template': "Duration Interaction ({:.1f}s) >> Duration Investigation ({:.1f}s)"
        },
        'A': {
            'description': 'Inefficient Feeding: 4-3 >> 6-5. Long interactions but little feeding.',
            'cycle_time_s': 155, # same long visit time
            'prob_interact': 0.9, # Not picky
            'prob_feed': 0.1,     # Low likelihood of feeding during long interaction
            'equation_template': "Duration Interaction ({:.1f}s) >> Duration Feeding ({:.1f}s)"
        }
    }
    
    results = {}

    for key, params in scenarios.items():
        print(f"--- Simulating Scenario for Answer Choice {key} ---")
        print(params['description'])
        
        # An insect can only start a new investigation after its previous cycle (visit + travel) is over.
        num_investigations = math.floor(total_time_seconds / params['cycle_time_s'])
        
        num_interactions = math.floor(num_investigations * params['prob_interact'])
        num_feedings = math.floor(num_interactions * params['prob_feed'])
        
        # Ethogram event counts (n) per hour
        n1_per_hour = num_investigations
        n3_per_hour = num_interactions
        n5_per_hour = num_feedings

        fitness_score = (num_interactions * 1) + (num_feedings * 5)
        results[key] = fitness_score
        
        print(f"Calculated events in 1 hour: n(1) Investigations={n1_per_hour}, n(3) Interactions={n3_per_hour}, n(5) Feedings={n5_per_hour}")

        # Display the equation with calculated numbers
        if key in ['F', 'E']:
            print(f"Final Equation: {params['equation_template'].format(n3_per_hour, n1_per_hour)}")
        elif key == 'C':
            # dur(interaction) = cycle_time - investigation_time - travel_time = 155 - 5 - 30 = 120s
            # dur(investigation) = 5s
            print(f"Final Equation: {params['equation_template'].format(120.0, 5.0)}")
        elif key == 'A':
            # dur(interaction) = 120s
            # dur(feeding) = dur(interaction) * prob_feed (as a proxy) = 120 * 0.1 = 12s
            print(f"Final Equation: {params['equation_template'].format(120.0, 12.0)}")

        print(f"Resulting Plant Fitness Score: {fitness_score}\n")

    best_choice = max(results, key=results.get)
    print("--- Conclusion ---")
    print(f"The scenario with the highest fitness score is Choice {best_choice} with a score of {results[best_choice]}.")
    print("This pattern represents a high number of efficient visits, which is most beneficial for cross-pollination.")

# Execute the simulation
simulate_behavior_and_calculate_fitness()
<<<F>>>
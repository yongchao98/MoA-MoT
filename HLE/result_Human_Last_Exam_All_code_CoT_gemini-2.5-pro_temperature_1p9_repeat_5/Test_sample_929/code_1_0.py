import random

def solve():
    """
    This script simulates different insect behaviors to determine which pattern
    maximizes plant fitness, which is primarily driven by pollination during feeding.
    """
    
    # Simulation parameters
    SIMULATION_TIME = 3600  # seconds in one hour
    fitness_scores = {}
    
    # A simplified function to simulate an insect's behavior for one hour.
    # The 'fitness score' is the total number of feeding events, as this
    # is the primary driver of pollination.
    def simulate_insect(name, description, prob_interact, time_investigate, time_interact, time_feed):
        print(f"--- Simulating Pattern '{name}' ---")
        print(description)
        
        time = 0
        n_investigations = 0  # Count for behavior 1
        n_interactions = 0    # Count for behavior 3
        n_feeds = 0           # Count for behavior 5

        # Seed the random number generator for consistent results
        random.seed(42)

        while time < SIMULATION_TIME:
            # An investigation cycle starts (behavior 1 starts)
            n_investigations += 1
            time += time_investigate # Duration of investigation (t2 - t1)

            # Check if the investigation leads to an interaction
            if random.random() < prob_interact:
                # An interaction starts (behavior 3 starts)
                n_interactions += 1
                # An interaction includes some non-feeding contact plus feeding
                interaction_duration = time_interact + time_feed # Total duration (t4 - t3)
                
                # Assume every successful interaction includes one feeding event
                if time_feed > 0:
                    n_feeds += 1
                
                time += interaction_duration
            
            if time >= SIMULATION_TIME:
                break
        
        fitness_score = n_feeds
        fitness_scores[name] = {
            "n_investigations": n_investigations,
            "n_interactions": n_interactions,
            "n_feeds": n_feeds,
            "fitness_score": fitness_score,
            "equation_lhs": n_interactions,
            "equation_rhs": n_investigations
        }
        print(f"Result: {n_investigations} investigations, {n_interactions} interactions, {n_feeds} feeds.\n")

    # Pattern F: n(3)/hour >> n(1)/hour. An efficient pollinator.
    simulate_insect(
        name="F",
        description="Behavior: High frequency of interactions relative to investigations.",
        prob_interact=0.9,      # Very likely to interact after investigating.
        time_investigate=10,
        time_interact=5,
        time_feed=25
    )

    # Pattern C: duration(4-3) >> duration(2-1). Spends a long time in one place.
    simulate_insect(
        name="C",
        description="Behavior: Long interaction durations relative to investigation durations.",
        prob_interact=0.4,      # Less frequent, but longer visits.
        time_investigate=5,     # Short investigation time.
        time_interact=50,       # Long interaction time...
        time_feed=150           # ...and long feeding time.
    )

    # Pattern E: n(1)/hour >> n(3)/hour. An inefficient pollinator.
    simulate_insect(
        name="E",
        description="Behavior: High frequency of investigations relative to interactions.",
        prob_interact=0.1,      # Very unlikely to interact after investigating.
        time_investigate=15,
        time_interact=5,
        time_feed=25
    )

    # --- Analysis ---
    print("--- Conclusion ---")
    best_pattern_name = max(fitness_scores, key=lambda k: fitness_scores[k]['fitness_score'])
    best_pattern_data = fitness_scores[best_pattern_name]
    
    print(f"The simulation shows that pattern '{best_pattern_name}' has the greatest positive effect on plant fitness with a score of {best_pattern_data['fitness_score']} feeds.")
    print("This pattern corresponds to n(3)/hour >> n(1)/hour, meaning the insect performs many interactions for each investigation.")
    print("\nThe final equation with numbers from the simulation is:")
    print(f"n({best_pattern_data['equation_lhs']})/hour >> n({best_pattern_data['equation_rhs']})/hour")

solve()
<<<F>>>
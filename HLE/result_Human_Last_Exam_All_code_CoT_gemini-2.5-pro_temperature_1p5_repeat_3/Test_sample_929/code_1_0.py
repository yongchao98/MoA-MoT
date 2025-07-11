import math

def solve_pollination_problem():
    """
    Analyzes different insect behavioral patterns to find which one most
    positively affects plant fitness, measured by pollinated flowers.
    """
    
    # --- Model Constants ---
    # The number of flowers on the plant's umbel available for pollination.
    NUM_FLOWERS_ON_PLANT = 30
    # The probability that one second of feeding results in pollination.
    POLLINATION_PROB_PER_SECOND = 0.25

    def calculate_plant_fitness(scenario_letter, description, n_investigations, n_interactions, dur_interaction, dur_feeding):
        """
        Calculates a fitness score based on insect behavior.
        
        Fitness = (Number of Unique Flowers Visited) * (Chance of Pollination per Visit)
        """
        print(f"--- Scenario {scenario_letter}: {description} ---")
        
        # --- Rule Checks ---
        # A feeding event (5-6) is a type of interaction (3-4).
        # Therefore, n(feeding) cannot be > n(interactions) and dur(feeding) cannot be > dur(interaction).
        if dur_feeding > dur_interaction:
            print("Logical impossibility: Feeding duration cannot exceed interaction duration.")
            print("Fitness Score: 0\n")
            return 0
        # This check is for option D, though it's covered by how we structure the simulation.
        # if n_feeding > n_interactions, it would also be impossible.

        # --- Calculations ---
        # 1. Number of Unique Flowers Visited
        # Assumes each interaction is an attempt on a new, random flower.
        # The number of unique flowers visited is capped by the total number of interactions
        # and the number of flowers available on the plant.
        flowers_visited = min(n_interactions, NUM_FLOWERS_ON_PLANT)

        # 2. Chance of Pollination per Visit
        # The probability of pollinating a single flower during a feeding bout.
        # P(success) = 1 - P(failure)
        # P(failure) = P(no pollination in 1s) ^ duration
        chance_of_pollination = 1 - (1 - POLLINATION_PROB_PER_SECOND)**dur_feeding

        # 3. Final Fitness Score
        # The total expected number of flowers pollinated in one hour.
        fitness_score = flowers_visited * chance_of_pollination

        print(f"Equation: Fitness = (Flowers Visited) * (Chance of Pollination per Visit)")
        print(f"Numbers:  Fitness = {flowers_visited} * {chance_of_pollination:.2f}")
        print(f"Fitness Score: {fitness_score:.2f}\n")
        return fitness_score

    # Define baseline parameters for a "typical" insect visit
    base_n_invest = 10
    base_n_interact = 10
    base_dur_interact = 10
    base_dur_feed = 8

    # --- Evaluate Each Answer Choice ---

    # A. Interaction duration is much greater than feeding duration.
    #    (Lots of non-feeding contact, inefficient for pollination)
    calculate_plant_fitness(
        "A", "4-3 >> 6-5",
        n_investigations=base_n_invest, n_interactions=base_n_interact,
        dur_interaction=50, dur_feeding=2)

    # B. Feeding duration is much greater than interaction duration.
    #    (Logically impossible)
    calculate_plant_fitness(
        "B", "6-5 >> 4-3",
        n_investigations=base_n_invest, n_interactions=base_n_interact,
        dur_interaction=10, dur_feeding=50)

    # C. Interaction duration is much greater than investigation duration.
    #    (Long contact time, but not necessarily many different flowers)
    calculate_plant_fitness(
        "C", "4-3 >> 2-1",
        n_investigations=base_n_invest, n_interactions=base_n_interact,
        dur_interaction=50, dur_feeding=40)

    # D. Number of feeding events >> number of interaction events.
    #    (Logically impossible, not directly simulated but would yield 0)
    print("--- Scenario D: n(5)/hour >> n(3)/hour ---")
    print("Logical impossibility: The number of feeding events cannot exceed the number of interactions.")
    print("Fitness Score: 0\n")
    
    # E. Number of investigations >> number of interactions.
    #    (Many approaches, but very few contacts. Poor for pollination)
    calculate_plant_fitness(
        "E", "n(1)/hour >> n(3)/hour",
        n_investigations=50, n_interactions=5,
        dur_interaction=base_dur_interact, dur_feeding=base_dur_feed)

    # F. Number of interactions >> number of investigations.
    #    (Insect efficiently visits many flowers on each approach. Excellent for pollination)
    calculate_plant_fitness(
        "F", "n(3)/hour >> n(1)/hour",
        n_investigations=5, n_interactions=50, # Note n_interactions > NUM_FLOWERS
        dur_interaction=base_dur_interact, dur_feeding=base_dur_feed)

if __name__ == '__main__':
    solve_pollination_problem()
import sys

def simulate_bacterial_fitness():
    """
    Simulates bacterial population dynamics to explain rapid resistance acquisition.

    This model demonstrates the principle from answer choice B. A rare resistance
    mutation often has a fitness cost, preventing it from spreading. However,
    a subsequent compensatory mutation can drastically increase fitness,
    allowing the new strain to rapidly take over the population.
    """

    # --- Step 1: Define fitness values ---
    # These are the key numbers in our population growth equation.
    # Fitness is a measure of reproductive success in the given environment.
    fitness_wt = 1.0   # Wild Type: Baseline fitness
    fitness_r = 0.8    # Resistance-only Mutant (R): Lower fitness due to cost
    fitness_cr = 1.5   # Compensatory + Resistance Mutant (CR): Very high fitness

    # --- Step 2: Initialize the population ---
    # Start with a large population of Wild Type bacteria.
    # Using dictionaries to keep track of counts.
    population = {'WT': 1000000, 'R': 0, 'CR': 0}
    generations = 30

    # --- Step 3: Print the model parameters ---
    print("--- Bacterial Population Simulation ---")
    print("This simulation shows how a highly fit mutant can rapidly dominate a population.")
    print("\nThe population growth is based on this equation: New_Count = Old_Count * Fitness")
    print("Here are the numbers used in the final equation for fitness:")
    print(f"Fitness of Wild Type (WT) = {fitness_wt}")
    print(f"Fitness of Resistance-only Mutant (R) = {fitness_r}")
    print(f"Fitness of Compensatory+Resistance Mutant (CR) = {fitness_cr}\n")
    print(f"Generation 0: 100% WT")

    # --- Step 4: Run the simulation over generations ---
    for gen in range(1, generations + 1):
        # Introduce a rare resistance mutation at generation 2
        if gen == 2:
            population['WT'] -= 1
            population['R'] += 1

        # Introduce a rare compensatory mutation in an R strain at generation 5
        if gen == 5:
            if population['R'] > 0:
                population['R'] -= 1
                population['CR'] += 1

        # Calculate the weighted count for each strain based on its fitness
        weighted_wt = population['WT'] * fitness_wt
        weighted_r = population['R'] * fitness_r
        weighted_cr = population['CR'] * fitness_cr

        # The sum of weighted counts represents the total reproductive output
        total_fitness_output = weighted_wt + weighted_r + weighted_cr
        if total_fitness_output == 0:
            break

        # The new proportion of each strain is its share of the total reproductive output
        total_population = sum(population.values())
        population['WT'] = round(total_population * (weighted_wt / total_fitness_output))
        population['R'] = round(total_population * (weighted_r / total_fitness_output))
        population['CR'] = round(total_population * (weighted_cr / total_fitness_output))

        # Print the population state at key moments to show the change
        if gen in [2, 5, 10, 15, 20, 30]:
            total_current_pop = sum(population.values())
            if total_current_pop > 0:
                percent_wt = round(100 * population['WT'] / total_current_pop, 2)
                percent_r = round(100 * population['R'] / total_current_pop, 2)
                percent_cr = round(100 * population['CR'] / total_current_pop, 2)
                print(f"Generation {gen:2}: WT: {percent_wt:>5}%, R: {percent_r:>5}%, CR: {percent_cr:>5}%")


    print("\n--- Conclusion ---")
    print("The R mutant, with its fitness cost, struggles to expand.")
    print("However, once the highly-fit CR mutant appears, it rapidly outcompetes all others.")
    print("This explains how a bacterium without lateral transfer can still acquire resistance at a fast pace.")

simulate_bacterial_fitness()
<<<B>>>
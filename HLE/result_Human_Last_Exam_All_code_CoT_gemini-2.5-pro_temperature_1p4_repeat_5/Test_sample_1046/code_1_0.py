def model_bacterial_fitness():
    """
    This script conceptually models the fitness changes in a bacterial population
    that acquires resistance through mutation, as described in the problem.
    """

    # --- Initial State ---
    # Let's assume a baseline fitness of 1.0 for the original, non-resistant bacterium.
    wild_type_fitness = 1.0
    print("Step 1: Start with a population of wild-type bacteria.")
    print(f"Initial fitness of a single wild-type bacterium = {wild_type_fitness}")
    print("-" * 60)

    # --- A Rare Resistance Mutation Occurs ---
    # This mutation often comes with a cost, reducing the bacterium's fitness.
    resistance_mutation_cost = 0.15
    resistant_mutant_fitness = wild_type_fitness - resistance_mutation_cost
    print("Step 2: A rare mutation occurs, granting drug resistance.")
    print("However, this mutation has a fitness cost.")
    print(f"Fitness of the new resistant mutant = {wild_type_fitness} - {resistance_mutation_cost} = {resistant_mutant_fitness}")
    print("This strain is less fit and would likely be outcompeted without the drug.")
    print("-" * 60)

    # --- A Compensatory Mutation Follows ---
    # A second mutation can compensate for the initial cost, restoring or even improving fitness.
    # This is the key insight from Option B.
    compensatory_mutation_gain = 0.25
    final_fitness = resistant_mutant_fitness + compensatory_mutation_gain
    print("Step 3: A second, compensatory mutation occurs in the resistant strain.")
    print("This new mutation increases fitness, making the strain highly competitive.")
    print(f"Final fitness of the compensated resistant mutant = {resistant_mutant_fitness} + {compensatory_mutation_gain} = {final_fitness:.2f}")
    print("-" * 60)

    # --- Conclusion ---
    print("Conclusion:")
    print(f"The final strain has a fitness of {final_fitness:.2f}, which is greater than the original wild-type's fitness of {wild_type_fitness}.")
    print("This high fitness allows it to rapidly spread through the population.")
    print("When you add the benefit of cross-resistance (one mutation resisting multiple drugs),")
    print("this mechanism can plausibly match the speed of lateral gene transfer.")
    print("This supports Option B as the best explanation.")

# Execute the function
model_bacterial_fitness()
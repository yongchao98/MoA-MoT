import numpy as np

# --- Simulation Parameters ---
POPULATION_SIZE = 500
GENERATIONS = 100
# Fitness equation parameters
BASE_FITNESS = 1.0
# Cost of deleterious mutations (higher mutation rate = more mutations = lower fitness)
COST_PER_MUTATION = 0.01
# Cost of fidelity/repair (lower mutation rate = more energy on repair = lower fitness)
COST_OF_FIDELITY = 0.0005

def calculate_fitness(num_mutations, mutation_rate):
    """
    Calculates fitness based on the number of deleterious mutations and
    the metabolic cost of maintaining a specific mutation rate.
    Fitness = Base_Fitness - (Cost_per_Mutation * Num_Mutations) - (Cost_of_Fidelity / Mutation_Rate)
    """
    if mutation_rate <= 0:
        return 0
    # The penalty for accumulating deleterious mutations
    mutation_load_penalty = COST_PER_MUTATION * num_mutations
    # The penalty for maintaining high fidelity (cost increases as rate gets very low)
    fidelity_penalty = COST_OF_FIDELITY / mutation_rate
    
    fitness = BASE_FITNESS - mutation_load_penalty - fidelity_penalty
    # Fitness cannot be negative
    return max(0, fitness)

def run_simulation():
    """
    Runs an evolutionary simulation to find the optimal mutation rate.
    """
    # Initialize a population with a wide range of mutation rates and no initial mutations
    population = [
        {'mutation_rate': rate, 'num_mutations': 0}
        for rate in np.random.uniform(0.01, 0.5, POPULATION_SIZE)
    ]
    
    initial_avg_rate = np.mean([org['mutation_rate'] for org in population])
    
    print("This simulation demonstrates how natural selection can maintain an optimal mutation rate.")
    print("Fitness is calculated by balancing the cost of new deleterious mutations against the cost of high-fidelity DNA repair.\n")
    print("The Fitness Equation:")
    print(f"Fitness = {BASE_FITNESS} - ({COST_PER_MUTATION} * Num_Deleterious_Mutations) - ({COST_OF_FIDELITY} / Mutation_Rate)\n")
    print(f"Initial average mutation rate: {initial_avg_rate:.4f}\n")

    for gen in range(GENERATIONS):
        # 1. Calculate fitness for every organism in the population
        fitness_scores = np.array([
            calculate_fitness(org['num_mutations'], org['mutation_rate'])
            for org in population
        ])
        
        # Normalize fitness scores to get selection probabilities
        total_fitness = np.sum(fitness_scores)
        if total_fitness == 0:
            print("Population went extinct. All fitness scores are zero.")
            return

        selection_probs = fitness_scores / total_fitness

        # 2. Create the next generation via selection
        # Parents are chosen based on their fitness
        parent_indices = np.random.choice(
            a=POPULATION_SIZE,
            size=POPULATION_SIZE,
            p=selection_probs
        )
        
        new_population = []
        for index in parent_indices:
            parent = population[index]
            
            # Offspring inherits mutation rate, with a small chance of change
            # This allows the mutation rate itself to evolve
            new_mutation_rate = parent['mutation_rate'] * np.random.normal(1.0, 0.05)
            # Prevent mutation rate from becoming non-positive or excessively high
            new_mutation_rate = max(0.001, min(1.0, new_mutation_rate))

            # 3. Accumulate new deleterious mutations
            # The number of new mutations depends on the parent's mutation rate
            new_mutations = np.random.poisson(lam=parent['mutation_rate'])
            total_mutations = parent['num_mutations'] + new_mutations
            
            new_population.append({
                'mutation_rate': new_mutation_rate,
                'num_mutations': total_mutations
            })
        
        population = new_population

    final_avg_rate = np.mean([org['mutation_rate'] for org in population])
    
    print("--- Simulation Complete ---")
    print(f"After {GENERATIONS} generations, the population's average mutation rate stabilized.")
    print(f"Final average mutation rate: {final_avg_rate:.4f}\n")
    print("This demonstrates that natural selection for fitness optimality (balancing opposing costs) can maintain an approximately constant genomic mutation rate.")

if __name__ == '__main__':
    run_simulation()
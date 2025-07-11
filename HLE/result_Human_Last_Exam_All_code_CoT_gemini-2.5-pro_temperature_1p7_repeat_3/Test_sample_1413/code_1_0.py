import numpy as np
from scipy.signal import convolve2d
import random

# Set a seed for reproducibility of the random process
random.seed(42)
np.random.seed(42)

# --- Configuration Constants ---

# Genetic Algorithm parameters
POPULATION_SIZE = 100
EVO_GENERATIONS = 50  # More generations can yield better results but take longer
ELITISM_RATE = 0.1    # Keep the top 10% of the population for the next generation
MUTATION_RATE = 0.015 # Chance for each cell to flip during mutation

# Game of Life simulation parameters
START_AREA_SIZE = 12
# The simulation grid must be larger to contain pattern expansion
SIM_GRID_SIZE = 150
# Max simulation steps to check for stability. Some patterns take thousands of steps.
MAX_GOL_STEPS = 6000
TARGET_FINAL_POPULATION = 100

def run_simulation(initial_pattern):
    """
    Runs a Conway's Game of Life simulation for a given initial pattern.
    Returns the final population count if it stabilizes above the target, otherwise 0.
    """
    # Place the 12x12 pattern in the center of a larger simulation grid
    grid = np.zeros((SIM_GRID_SIZE, SIM_GRID_SIZE), dtype=np.int8)
    start_pos = (SIM_GRID_SIZE - START_AREA_SIZE) // 2
    grid[start_pos:start_pos + START_AREA_SIZE, start_pos:start_pos + START_AREA_SIZE] = initial_pattern

    history = {grid.tobytes()}
    kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])

    for step in range(MAX_GOL_STEPS):
        # Use convolution to count live neighbors for each cell efficiently
        # 'wrap' mode simulates an infinite grid
        num_neighbors = convolve2d(grid, kernel, mode='same', boundary='wrap')
        
        # Apply Game of Life rules
        born = (num_neighbors == 3) & (grid == 0)
        survive = ((num_neighbors == 2) | (num_neighbors == 3)) & (grid == 1)
        
        new_grid = np.zeros_like(grid)
        new_grid[born | survive] = 1

        # Check for stability
        if np.array_equal(grid, new_grid): # Still life
            return np.sum(new_grid)
            
        new_grid_hash = new_grid.tobytes()
        if new_grid_hash in history: # Oscillator / repeating pattern
            return np.sum(new_grid)
        
        history.add(new_grid_hash)
        grid = new_grid
        
        if np.sum(grid) == 0: # Died out
            return 0
            
    return 0 # Did not stabilize within max steps

def calculate_fitness(individual):
    """
    Calculates the fitness of an individual. The fitness is the initial
    number of cells if the final stable population is > 100, otherwise 0.
    """
    final_pop = run_simulation(individual)
    if final_pop > TARGET_FINAL_POPULATION:
        return np.sum(individual)
    return 0

def create_initial_population():
    """Creates the starting population, including random and known patterns."""
    population = []

    # 1. Add known high-performing methuselahs to seed the search
    # Acorn (7 cells -> 633)
    acorn = np.zeros((START_AREA_SIZE, START_AREA_SIZE), dtype=np.int8)
    acorn[5, 4] = 1; acorn[6, 6] = 1; acorn[7, 3:5] = 1; acorn[7, 6:9] = 1
    population.append(acorn)

    # R-pentomino (5 cells -> 116)
    r_pento = np.zeros((START_AREA_SIZE, START_AREA_SIZE), dtype=np.int8)
    r_pento[5, 5:7] = 1; r_pento[6, 4:6] = 1; r_pento[7, 5] = 1
    population.append(r_pento)

    # 2. Fill the rest with random individuals
    while len(population) < POPULATION_SIZE:
        # Create random patterns with a 20-50% density
        density = random.uniform(0.2, 0.5)
        individual = np.random.choice([0, 1], size=(START_AREA_SIZE, START_AREA_SIZE), p=[1 - density, density])
        population.append(individual.astype(np.int8))
        
    return population

def select_parent(sorted_population):
    """Selects a parent using tournament selection."""
    tournament_size = 5
    participants = random.choices(sorted_population, k=tournament_size)
    # The list is already sorted, so the first element is the best
    return participants[0][1]

def crossover(parent1, parent2):
    """Performs single-point crossover on two parents."""
    p1_flat = parent1.flatten()
    p2_flat = parent2.flatten()
    crossover_point = random.randint(1, len(p1_flat) - 2)
    child_flat = np.concatenate((p1_flat[:crossover_point], p2_flat[crossover_point:]))
    return child_flat.reshape((START_AREA_SIZE, START_AREA_SIZE))

def mutate(individual):
    """Applies mutation to an individual by flipping random cells."""
    mask = np.random.rand(*individual.shape) < MUTATION_RATE
    individual[mask] = 1 - individual[mask]
    return individual

def main():
    """Main function to run the genetic algorithm."""
    population = create_initial_population()
    best_overall_fitness = 0
    best_pattern = None

    print(f"Starting genetic algorithm search...")
    print(f"Population size: {POPULATION_SIZE}, Generations: {EVO_GENERATIONS}")
    print("-" * 30)

    for gen in range(EVO_GENERATIONS):
        fitness_scores = [calculate_fitness(ind) for ind in population]
        
        # Pair individuals with their fitness and sort by fitness (descending)
        pop_with_fitness = sorted(zip(fitness_scores, population), key=lambda x: x[0], reverse=True)

        current_best_fitness = pop_with_fitness[0][0]
        if current_best_fitness > best_overall_fitness:
            best_overall_fitness = current_best_fitness
            best_pattern = pop_with_fitness[0][1]
            # Since fitness is initial population, this is our target value
            final_pop_for_best = run_simulation(best_pattern)
            print(f"Gen {gen+1}: New best found! Initial Cells: {best_overall_fitness}, Stabilizes to: {final_pop_for_best} cells")

        # Create the next generation
        new_population = []
        
        # Elitism: carry over the best individuals
        elite_count = int(POPULATION_SIZE * ELITISM_RATE)
        elites = [ind for fitness, ind in pop_with_fitness[:elite_count]]
        new_population.extend(elites)

        # Create the rest of the new population through crossover and mutation
        while len(new_population) < POPULATION_SIZE:
            parent1 = select_parent(pop_with_fitness)
            parent2 = select_parent(pop_with_fitness)
            child = crossover(parent1, parent2)
            child = mutate(child)
            new_population.append(child)
            
        population = new_population

    print("-" * 30)
    print("Search complete.")

    if best_pattern is not None:
        initial_cells = np.sum(best_pattern)
        final_cells = run_simulation(best_pattern)
        print(f"The best pattern found has an initial population of {initial_cells} cells.")
        print(f"It evolves for thousands of steps and eventually stabilizes at {final_cells} cells.")
        print("\nFinal Answer (Greatest number of initial live cells found):")
        print(initial_cells)
        # Final answer format for the platform
        print(f"<<<{initial_cells}>>>")

    else:
        print("Could not find a pattern that stabilizes to over 100 cells.")
        print("Try increasing POPULATION_SIZE or EVO_GENERATIONS.")
        print(f"<<<0>>>")

if __name__ == '__main__':
    main()
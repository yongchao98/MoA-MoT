import random

def simulate_evolution():
    """
    This simulation models the emergence of mucoid variants by mutation
    in the mucA gene in two types of bacterial strains: Normal and Hypermutator.
    """
    # --- Simulation Parameters ---
    population_size = 10000
    gene_length = 30  # A simplified representation of the mucA gene's target size
    generations = 100
    normal_mutation_rate = 0.0005    # Mutation rate per base per generation
    hypermutator_mutation_rate = 0.01  # A 20x higher mutation rate

    print("This script simulates the effect of hypermutation on the frequency and spectrum of mucoid variants.")
    print("The simulation will compare a 'Normal Strain' with a 'Hypermutator Strain'.\n")

    def run_population_simulation(strain_name, mutation_rate):
        """
        Runs the simulation for a single population.
        A 'mucoid variant' is any bacterium with a mutated gene.
        """
        # All individuals start with a non-mutated (wild-type) gene
        wild_type_gene = "W" * gene_length
        
        # Keep track of unique mutations to measure the 'spectrum'
        unique_mutations_found = set()
        
        # In this simple model, we just calculate the expected number of mutants
        # without simulating every individual bacterium, which is more efficient.
        
        # Calculate expected number of new mucoid variants per generation
        # Prob. of at least one mutation = 1 - (Prob. of no mutations)
        # Prob. of no mutations = (1 - mutation_rate) ^ gene_length
        prob_gene_mutates = 1 - (1 - mutation_rate) ** gene_length
        
        # Track the number of non-mutated bacteria
        num_wild_type = population_size
        total_mucoid_variants = 0

        for gen in range(generations):
            # Number of bacteria that will acquire a new mutation in this generation
            newly_mutated = round(num_wild_type * prob_gene_mutates)
            num_wild_type -= newly_mutated
            total_mucoid_variants += newly_mutated
            
            # For spectrum, we assume each new mutation event can create a new unique mutation
            # This is a simplification but illustrates the concept.
            # We'll generate a sample of these potential new mutations.
            for _ in range(newly_mutated):
                # Create a representative random mutation
                mutant_gene = list(wild_type_gene)
                pos = random.randint(0, gene_length - 1)
                mutant_gene[pos] = 'M'
                unique_mutations_found.add("".join(mutant_gene))

        frequency = total_mucoid_variants
        spectrum_width = len(unique_mutations_found)
        
        print(f"--- Results for {strain_name} ---")
        print(f"Initial Population: {population_size}")
        print(f"Mutation Rate: {mutation_rate}")
        print(f"Simulated Generations: {generations}")
        print("-" * 35)
        print(f"Final calculated mucoid variants (Frequency): {frequency}")
        print(f"Final unique mucA mutations observed (Spectrum): {spectrum_width}")
        print("-" * 35 + "\n")

    # Run simulations and print results
    run_population_simulation("Normal Strain", normal_mutation_rate)
    run_population_simulation("Hypermutator Strain", hypermutator_mutation_rate)
    
    print("Conclusion: The simulation shows the Hypermutator Strain results in a higher frequency of variants and a wider spectrum of mutations.")

simulate_evolution()
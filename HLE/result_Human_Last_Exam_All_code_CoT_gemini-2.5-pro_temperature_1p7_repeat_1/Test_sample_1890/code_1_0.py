import random

def run_bacterial_evolution_simulation():
    """
    Simulates the emergence of mucoid variants in Pseudomonas to answer the question.
    This model demonstrates the effects of different mutation rates on the frequency
    and variety of mutations in a target gene (mucA).
    """
    
    # --- Parameters ---
    # Our model gene. 'F' = functional. Any other character = mutated/non-functional.
    wild_type_mucA_gene = 'FFFFFFFFFF' 
    gene_length = len(wild_type_mucA_gene)
    
    # Simulation settings
    population_size = 5000
    generations = 100
    
    # Mutation rates per base per generation
    wild_type_rate = 0.0001
    hypermutator_rate = 0.005 # 50x higher mutation rate
    
    print("Plan: Simulate two bacterial populations (Wild-Type and Hypermutator).")
    print("We will track how many mucoid variants arise and how many different types of mutations cause them.\n")

    # --- Simulation Function ---
    def simulate_strain(strain_name, mutation_rate):
        print(f"--- Simulating: {strain_name} Strain ---")
        mucoid_variants_found = 0
        # A set stores only unique items, perfect for tracking the spectrum of mutations
        mutation_spectrum = set()
        
        for gen in range(generations):
            for _ in range(population_size):
                # Assume each bacterium starts with a functional mucA gene
                current_gene = list(wild_type_mucA_gene)
                mutated = False
                for i in range(gene_length):
                    if random.random() < mutation_rate:
                        current_gene[i] = 'M' # 'M' for mutated
                        mutated = True
                
                if mutated:
                    mucoid_variants_found += 1
                    # Record the specific mutation pattern
                    mutation_spectrum.add("".join(current_gene))

        total_simulated_bacteria = population_size * generations
        frequency = mucoid_variants_found / total_simulated_bacteria
        spectrum_width = len(mutation_spectrum)
        
        print(f"Result for {strain_name}:")
        print(f"Total Mucoid Variants Found: {mucoid_variants_found}")
        # Final equation with numbers as requested
        print(f"Frequency Calculation: {mucoid_variants_found} / ({population_size} * {generations}) = {frequency:.6f}")
        print(f"Spectrum of Mutations (variety): {spectrum_width}\n")
        
        return frequency, spectrum_width

    # --- Run Simulations and Conclude ---
    wt_freq, wt_spec = simulate_strain("Wild-Type", wild_type_rate)
    hm_freq, hm_spec = simulate_strain("Hypermutator", hypermutator_rate)

    print("--- Conclusion ---")
    if hm_freq > wt_freq:
        print(f"The frequency of mucoid variants is higher in Hypermutators ({hm_freq:.4f}) than in Wild-Type ({wt_freq:.4f}). -> INCREASED FREQUENCY")
    
    if hm_spec > wt_spec:
        print(f"The spectrum of mutations is wider in Hypermutators ({hm_spec}) than in Wild-Type ({wt_spec}). -> WIDER SPECTRUM")

    print("\nTherefore, hypermutator strains provide an increased frequency of mucoid variants with a wider spectrum of mucA mutations.")

# Run the simulation and print the results
run_bacterial_evolution_simulation()
<<<B>>>
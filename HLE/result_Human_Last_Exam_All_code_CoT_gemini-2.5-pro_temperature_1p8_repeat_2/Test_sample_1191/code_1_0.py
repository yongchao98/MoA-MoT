import math

def calculate_mutational_load():
    """
    This function models the effect of the genomic mutation rate (U) on
    the mean fitness of a population.
    
    The concept is that most new mutations are deleterious. A higher mutation
    rate means more deleterious mutations per individual, which lowers the
    overall population fitness. This creates a selective pressure to reduce
    the mutation rate.

    The model used is W = e^(-U), where:
    - W is the mean population fitness.
    - U is the genomic rate of deleterious mutations.
    - 'e' is the base of the natural logarithm.
    """
    
    print("Modeling the relationship between genomic mutation rate (U) and population fitness (W).")
    print("The equation is: W = e^(-U)")
    print("-" * 50)
    
    # A range of possible deleterious mutation rates per genome per generation.
    genomic_mutation_rates = [0.01, 0.1, 0.5, 1.0, 2.0]
    
    for U in genomic_mutation_rates:
        # Calculate the mean fitness of the population.
        # As U increases, the fitness W decreases.
        mean_fitness = math.exp(-U)
        
        print(f"When the deleterious mutation rate (U) is {U:.2f}:")
        print(f"The equation becomes: W = e^(-{U:.2f})")
        print(f"The resulting mean population fitness (W) is {mean_fitness:.4f}")
        print("This lower fitness is what natural selection acts against, favoring a lower U.\n")

# Run the calculation and print the results
calculate_mutational_load()
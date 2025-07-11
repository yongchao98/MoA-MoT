import textwrap

def explain_mutation_rate_stability():
    """
    Explains the reasoning behind the dominant factor maintaining a constant
    genomic mutation rate and prints the final answer.
    """
    # Introduction to the problem
    print("Analyzing the factors that maintain a constant genomic mutation rate:")
    print("-" * 60)
    print("The key insight is that the mutation rate itself is a biological trait that is subject to evolutionary forces.\n")

    # Analysis of each option
    options = {
        'A': "Natural selection for fitness optimality: A very high mutation rate creates a 'mutational load' (too many harmful mutations), which is selected against. A very low mutation rate is metabolically expensive to maintain (high-fidelity DNA repair costs energy) and can limit adaptive potential. Therefore, natural selection provides a stabilizing pressure, favoring an intermediate rate and keeping it relatively constant.",
        'B': "Genetic drift in small populations: Genetic drift is a random force that causes allele frequencies to fluctuate. It would lead to random wandering of the mutation rate, not stability.",
        'C': "Equilibrium between beneficial and deleterious mutations: This equilibrium describes the consequence of a given mutation rate on a population's fitness, but it doesn't regulate the underlying rate at which mutations are generated.",
        'D': "Homogeneous mutation distribution across genomic sites: This describes the spatial pattern of where mutations occur, not their overall frequency or rate. This premise is also generally false, as mutation hotspots exist.",
        'E': "The stochastic nature of mutational events: 'Stochastic' simply means random. The randomness of individual mutational events does not, by itself, explain the stability of their overall rate over evolutionary time."
    }

    print("Evaluating the answer choices:\n")
    for option, explanation in options.items():
        print(f"Option {option}:")
        # Use textwrap for clean printing
        wrapped_explanation = textwrap.fill(explanation, width=70, initial_indent="  ", subsequent_indent="  ")
        print(wrapped_explanation)
        print()

    # Final conclusion
    print("-" * 60)
    print("Conclusion: The most robust explanation is that natural selection stabilizes the mutation rate by selecting against both extremes. This balances the cost of deleterious mutations with the cost of high-fidelity DNA maintenance.")
    
    final_answer = 'A'
    print(f"\nTherefore, the dominant factor is: {final_answer}")

# Run the explanation function
explain_mutation_rate_stability()
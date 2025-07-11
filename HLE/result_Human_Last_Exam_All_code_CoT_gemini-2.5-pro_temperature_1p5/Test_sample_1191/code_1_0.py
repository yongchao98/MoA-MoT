import textwrap

def explain_mutation_rate_equilibrium():
    """
    This function explains the evolutionary principle that maintains a constant genomic mutation rate.
    """
    # The core concept is that the mutation rate itself is an evolvable trait
    # shaped by competing selective pressures.

    # Pressure against a high mutation rate:
    cost_of_high_rate = "A high mutation rate introduces a large number of deleterious mutations, creating a 'genetic load' that reduces the average fitness of the population."

    # Pressure against a very low mutation rate:
    cost_of_low_rate = "A very low mutation rate limits the potential for adaptation by reducing the supply of beneficial mutations. It can also be metabolically expensive to achieve near-perfect replication fidelity."

    # The resolution of these competing pressures is an equilibrium.
    conclusion = "The observed constant genomic mutation rate is an evolutionary equilibrium, balancing the immediate harm from deleterious mutations with the long-term need for beneficial mutations."

    # Map this conclusion to the answer choices.
    answer_choice = 'C'
    answer_text = "Equilibrium between beneficial and deleterious mutations"

    print("Analyzing the forces that shape the genomic mutation rate:")
    print("-" * 60)
    print("Factor 1: The Cost of Deleterious Mutations")
    print(textwrap.fill(cost_of_high_rate, width=60))
    print("\nFactor 2: The Need for Beneficial Mutations & Cost of Fidelity")
    print(textwrap.fill(cost_of_low_rate, width=60))
    print("-" * 60)
    print("\nConclusion:")
    print(textwrap.fill(conclusion, width=60))
    print(f"\nThis corresponds to Answer Choice '{answer_choice}': {answer_text}")

explain_mutation_rate_equilibrium()
import textwrap

def explain_mutation_rate_equilibrium():
    """
    Explains the evolutionary reasoning for a stable genomic mutation rate.
    """
    print("The question asks which factor dominates to maintain a constant genomic mutation rate under mutation pressure.")
    print("-" * 80)
    
    explanation = """
    The maintenance of a specific, relatively constant genomic mutation rate across diverse life forms is a classic topic in evolutionary biology. It is best understood as an evolutionary trade-off, leading to an equilibrium.

    1.  The Problem with a HIGH Mutation Rate:
        Mutations are frequently neutral or slightly harmful (deleterious). A very high mutation rate would introduce a large number of deleterious mutations into the population in each generation. This increases the 'mutational load', reducing the average fitness of the population and potentially leading to extinction (a concept called 'mutational meltdown').

    2.  The Problem with a LOW Mutation Rate:
        While avoiding deleterious mutations is beneficial, a mutation rate that is too low is also disadvantageous. Beneficial mutations are the raw material for adaptation. If the rate is too low, beneficial mutations arise too infrequently, crippling a population's ability to adapt to new challenges, such as changing climates, new diseases, or competing species.

    3.  The Equilibrium Solution:
        Therefore, the mutation rate itself is subject to natural selection. The rate that tends to persist is an evolutionary compromise—an equilibrium point. This rate is high enough to generate sufficient beneficial mutations for adaptation but not so high that the burden of deleterious mutations becomes unsustainable. This balance is best described as an 'Equilibrium between beneficial and deleterious mutations'.
    """
    
    # Use textwrap to format the explanation nicely
    wrapped_explanation = textwrap.dedent(explanation).strip()
    print(wrapped_explanation)
    print("-" * 80)
    print("Conclusion: Option C best captures this fundamental evolutionary trade-off.")

# Execute the explanation function
explain_mutation_rate_equilibrium()
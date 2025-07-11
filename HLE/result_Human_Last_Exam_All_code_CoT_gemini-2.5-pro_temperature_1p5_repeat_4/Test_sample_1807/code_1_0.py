import textwrap

def explain_sex_linked_fst():
    """
    Explains why sex determination systems can cause high Fst between males and females.
    """
    print("### Analysis of Genetic Differentiation (Fst) Between Sexes ###\n")

    # Define the observed phenomenon
    print("Observation: High Fst is found for SOME genetic markers when comparing males and females of a single population.")
    print("This means that for these specific markers, allele frequencies differ significantly between the sexes.\n")

    print("--- Evaluating the Potential Explanations ---\n")

    # Dictionary of explanations for each choice
    explanations = {
        'A. Genetic load': "Incorrect. Genetic load is the reduced fitness of a population due to deleterious alleles. While some deleterious alleles might be on sex chromosomes, genetic load itself is a measure of fitness, not a direct cause of allele frequency differentiation (Fst).",
        'B. XY vs ZW sex determining systems': "Correct. This is the most direct cause. In these systems, sexes have different chromosome compositions:\n"
                                           "   - XY System (e.g., humans): Males are XY, Females are XX. Any marker on the Y chromosome will only be found in males, leading to an Fst of 1. Markers on the X chromosome also differ in copy number (1 in males, 2 in females), which can cause Fst > 0.\n"
                                           "   - ZW System (e.g., birds): Females are ZW, Males are ZZ. Any marker on the W chromosome will only be found in females, leading to an Fst of 1.\n"
                                           "   - In both cases, markers on autosomes (non-sex chromosomes) would show Fst ≈ 0. This perfectly explains why only 'some markers' exhibit high differentiation.",
        'C. Reproductive isolation': "Incorrect. Males and females within a single population are, by definition, NOT reproductively isolated. They must interbreed for the population to exist.",
        'D. Local adaptation to thermal environments': "Incorrect. Local adaptation explains differentiation between populations living in different environments. It does not explain genetic differentiation between sexes that share the same environment and interbreed.",
        'E. Hybrid zone dynamics': "Incorrect. A hybrid zone is where two genetically distinct populations meet and interbreed. The question describes a single, panmictic population."
    }

    for option, reason in explanations.items():
        # Use textwrap for neat printing of long explanation strings
        wrapped_reason = textwrap.fill(reason, width=80, initial_indent='   ', subsequent_indent='   ')
        print(f"Choice {option}:\n{wrapped_reason}\n")

    print("--- Conclusion ---\n")
    print("The presence of sex chromosomes, which are inherited differently and are sometimes unique to one sex, is the primary reason for high Fst values at specific loci when comparing males and females.")

# Run the explanation function
explain_sex_linked_fst()
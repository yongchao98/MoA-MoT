def explain_gene_flow_effects():
    """
    Explains the effects of gene flow in a hybrid zone on various population genetic metrics.
    """
    print("This script explains why 'High Fst between populations' is the correct answer.")
    print("--------------------------------------------------------------------------\n")
    print("Background Concepts:")
    print("1. Hybrid Zone: A region where two genetically distinct populations meet and interbreed.")
    print("2. Gene Flow: The transfer of genetic material (alleles) from one population to another.\n")
    
    print("Analysis of Answer Choices:")
    
    # Choice A
    print("\nA. High Fst between populations")
    print("   - Fst measures population differentiation. An Fst of 0 means populations are identical, while 1 means they are completely different.")
    print("   - Gene flow is a homogenizing force; it mixes alleles between populations, making them more similar.")
    print("   - Therefore, the very process of gene flow acts to *reduce* Fst.")
    print("   - While the parent populations that form the hybrid zone must have had high Fst, sustained gene flow across the zone is fundamentally incompatible with maintaining high Fst. High gene flow and high Fst are opposing states.")
    print("   - Verdict: This is directly counteracted by gene flow.")

    # Choice B
    print("\nB. High Dxy between populations")
    print("   - Dxy measures the absolute genetic divergence between two populations (e.g., average number of nucleotide differences).")
    print("   - It reflects the amount of time the two populations have been separated.")
    print("   - Gene flow does not erase this history of divergence. The parent populations can have a high Dxy from a long period of separation before coming into contact.")
    print("   - Verdict: This can certainly occur.")

    # Choice C
    print("\nC. High Fis within a population")
    print("   - Fis is the inbreeding coefficient and measures the deficit or excess of heterozygotes within a population.")
    print("   - In a hybrid zone, pooling individuals from two distinct parental gene pools creates a statistical deficit of heterozygotes (known as the Wahlund effect).")
    print("   - This results in a high positive Fis value.")
    print("   - Verdict: This is an expected consequence of hybridization.")

    # Choice D
    print("\nD. High u (mutation rate) within a population")
    print("   - 'u' (or mu, μ) is the rate at which new mutations arise.")
    print("   - This is a fundamental molecular property of an organism's biology.")
    print("   - Gene flow is a population-level process and has no direct mechanistic effect on the molecular mutation rate.")
    print("   - A population in a hybrid zone could have a high mutation rate for unrelated reasons (e.g., exposure to a mutagen).")
    print("   - Verdict: This is independent of gene flow and can occur.")

    # Choice E
    print("\nE. High Pi within a population")
    print("   - Pi (π), or nucleotide diversity, measures the amount of genetic variation within a population.")
    print("   - By mixing alleles from two different, divergent source populations, gene flow increases the total amount of genetic variation within the hybrid zone itself.")
    print("   - Verdict: This is an expected consequence of gene flow into a hybrid zone.")

    print("\n--------------------------------------------------------------------------")
    print("Conclusion:")
    print("Gene flow's primary effect is to make populations more similar, thus decreasing Fst. A state of high gene flow cannot coexist with a state of high Fst. Therefore, 'High Fst' is the condition that cannot occur as a result of effective gene flow.")

explain_gene_flow_effects()
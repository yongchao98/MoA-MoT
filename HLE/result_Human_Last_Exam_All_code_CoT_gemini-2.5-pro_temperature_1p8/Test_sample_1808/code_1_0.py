def analyze_gene_flow_effects():
    """
    Analyzes which population genetic scenario is incompatible with
    gene flow across a hybrid zone.
    """
    print("Analyzing the effect of gene flow across a hybrid zone...")
    print("-" * 60)

    # Explanation of concepts
    print("Key Concepts:")
    print("  - Gene Flow: The movement of alleles between populations.")
    print("  - Hybrid Zone: A region where distinct populations meet and interbreed.")
    print("  - The primary effect of gene flow is homogenization, making populations more genetically similar.\n")

    # Evaluation of each option
    print("Evaluating the Options:")

    # Option A: High Fst
    print("A. High Fst between populations:")
    print("   Fst measures population differentiation. A high Fst (close to 1.0) means")
    print("   populations are very distinct with little genetic overlap. Gene flow directly")
    print("   counteracts this by mixing alleles, thus *lowering* Fst. Therefore,")
    print("   the presence of significant gene flow is inconsistent with maintaining a high Fst.")
    print("   This is a likely candidate for something that CANNOT occur.\n")

    # Option B: High Dxy
    print("B. High Dxy between populations:")
    print("   Dxy measures the absolute genetic divergence, reflecting historical separation time.")
    print("   Two populations can be highly diverged from their past (high Dxy) and still")
    print("   experience recent gene flow in a hybrid zone. Gene flow does not erase")
    print("   historical divergence quickly. So, this CAN occur.\n")

    # Option C: High Fis
    print("C. High Fis within a population:")
    print("   Fis measures inbreeding (excess of homozygotes). In a hybrid zone, if")
    print("   individuals prefer to mate with their own type (assortative mating) or if")
    print("   there is a mix of two distinct populations (Wahlund effect), it can lead")
    print("   to a deficit of heterozygotes, resulting in a high Fis. So, this CAN occur.\n")

    # Option D: High u (mutation rate)
    print("D. High μ (mutation rate) within a population:")
    print("   Mutation rate is a fundamental biological property of an organism's DNA.")
    print("   Gene flow is the movement of existing alleles and has no direct impact on the")
    print("   rate at which new mutations arise. So, this CAN occur.\n")

    # Option E: High Pi
    print("E. High Pi within a population:")
    print("   Pi (π) measures nucleotide diversity. Gene flow from a different population")
    print("   introduces new alleles, increasing the genetic variation within the hybrid zone.")
    print("   This directly leads to a higher Pi. So, this CAN occur.\n")

    print("-" * 60)
    print("Conclusion: Gene flow works to decrease Fst. Therefore, the one thing that")
    print("cannot occur *as a result of or in the presence of* significant gene flow")
    print("is a high Fst between the participating populations.")


if __name__ == "__main__":
    analyze_gene_flow_effects()
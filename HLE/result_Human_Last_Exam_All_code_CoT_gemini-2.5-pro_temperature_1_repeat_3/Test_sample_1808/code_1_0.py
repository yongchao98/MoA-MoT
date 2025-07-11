def analyze_hybrid_zone_scenarios():
    """
    This script provides a step-by-step analysis of the potential genetic outcomes
    when gene flow occurs across a hybrid zone, identifying the one that cannot occur.
    """
    print("Analysis of Population Genetics in a Hybrid Zone")
    print("="*50)
    print("Scenario: Gene flow is occurring across a hybrid zone, which is an area where two genetically distinct populations meet and interbreed.\n")

    print("Let's analyze each potential outcome:\n")

    print("A. High Fst between populations")
    print("   - Fst measures genetic differentiation between populations. A hybrid zone forms because two populations are already differentiated.")
    print("   - Therefore, a high Fst between the source populations is an expected starting condition.")
    print("   - Verdict: This CAN occur.\n")

    print("B. High Dxy between populations")
    print("   - Dxy measures the absolute genetic divergence between populations. Like Fst, it reflects the historical separation of the populations.")
    print("   - A high Dxy is expected between the two distinct populations that form the hybrid zone.")
    print("   - Verdict: This CAN occur.\n")

    print("C. High Fis within a population")
    print("   - Fis measures the heterozygote deficit within a population. A high positive Fis indicates fewer heterozygotes than expected under random mating (Hardy-Weinberg Equilibrium).")
    print("   - When you sample from a hybrid zone, you are pooling two genetically distinct groups. This is a classic example of the Wahlund effect, which causes a statistical deficit of heterozygotes.")
    print("   - Verdict: This CAN occur and is highly expected.\n")

    print("D. High u (mutation rate) within a population")
    print("   - 'u' (or μ) is the mutation rate, the fundamental rate at which new alleles arise from changes in DNA sequence.")
    print("   - Gene flow is the process of moving *existing* alleles between populations. It does not create new ones.")
    print("   - The mutation rate is a biochemical property of an organism and is not influenced by the population-level process of gene flow.")
    print("   - Verdict: This CANNOT occur as a result of gene flow.\n")

    print("E. High Pi within a population")
    print("   - Pi (π) is nucleotide diversity, a measure of genetic variation within a single population.")
    print("   - Gene flow introduces new alleles from one population into another. In a hybrid zone, the mixing of two different gene pools directly increases the amount of genetic variation.")
    print("   - This results in a higher Pi within the hybrid zone compared to either parent population.")
    print("   - Verdict: This CAN occur and is highly expected.\n")

    print("="*50)
    print("Conclusion: Gene flow affects the frequencies of existing alleles, while mutation creates new ones. These are independent processes. Therefore, gene flow does not cause a high mutation rate.")

if __name__ == '__main__':
    analyze_hybrid_zone_scenarios()

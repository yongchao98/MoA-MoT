import textwrap

def analyze_gwas_scenario():
    """
    Analyzes a multiple-choice question about misleading associations in GWAS
    due to linkage disequilibrium (LD).
    """

    print("Analyzing the question: Which combination of SNP alleles is most likely to provide a misleading association with a complex trait due to linkage disequilibrium?\n")
    
    # Core concepts
    print("--- Key Concepts ---")
    print("1. Complex Trait: Influenced by multiple genes, often on different chromosomes (polygenic).")
    print("2. Linkage Disequilibrium (LD): The tendency for alleles at nearby locations on a chromosome to be inherited together as a block.")
    print("3. Misleading Association: A strong statistical signal that points to a non-causal SNP (because it's in LD with the true causal SNP) or overemphasizes the importance of one genomic region, obscuring the true polygenic nature of the trait.\n")

    # Analysis of options
    print("--- Analysis of Answer Choices ---")

    # Option A
    analysis_a = "A. Three SNPs tightly linked within a single LD block on one chromosome.\n   - If a causal variant lies within this block, all three SNPs will show a strong statistical association.\n   - This creates a powerful, concentrated signal at one locus.\n   - For a polygenic trait caused by many scattered loci, this strong signal can overshadow other weaker, true signals, misleadingly suggesting this single region is the main driver.\n   - This is highly likely to be misleading about the overall genetic architecture."
    print(textwrap.indent(analysis_a, '  '))

    # Option B
    analysis_b = "B. Two SNPs located at the extreme ends of an LD block.\n   - This would also show association due to LD, but the signal is often weaker at the edges of a block. Less potent than scenario A."
    print(textwrap.indent(analysis_b, '  '))

    # Option C
    analysis_c = "C. Two SNPs located on different chromosomes.\n   - SNPs on different chromosomes are not in LD. An association with both would correctly suggest two independent loci, which is the opposite of a misleading result due to LD."
    print(textwrap.indent(analysis_c, '  '))

    # Option D
    analysis_d = "D. A single SNP located centrally within an LD block.\n   - This is a standard GWAS finding. The association is with the entire block. It doesn't create the same kind of overwhelming, concentrated signal as option A."
    print(textwrap.indent(analysis_d, '  '))
    
    # Option E
    analysis_e = "E. Three Tag SNPs predicting all alleles in an inherited haplotype.\n   - This describes the general methodology of using tag SNPs to represent an LD block (haplotype). Option A is a specific, potent *result* of this methodology that leads to a misleading conclusion."
    print(textwrap.indent(analysis_e, '  '))

    # Conclusion
    print("\n--- Conclusion ---")
    conclusion_text = "Option A is the most likely scenario. The strong, redundant signal from three tightly linked SNPs in a single block can mask the smaller effects from other true causal loci scattered across the genome. This misrepresents the complex, polygenic nature of the trait by making one locus appear overwhelmingly important."
    print(conclusion_text)

if __name__ == "__main__":
    analyze_gwas_scenario()
    print("\n<<<A>>>")

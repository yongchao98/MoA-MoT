import pandas as pd

def simulate_ld_association():
    """
    Simulates and illustrates how Tag SNPs can produce a misleading association
    with a trait due to Linkage Disequilibrium (LD).
    """
    print("### Simulation: Misleading Association due to LD and Tag SNPs ###\n")

    # Step 1: Define two haplotypes in an LD block. A haplotype is a sequence of
    # alleles on a single chromosome that are inherited together.
    # We will define one SNP as the true "Causal SNP" for a hypothetical disease.
    # Haplotype 1 is the 'Risk' haplotype and Haplotype 2 is the 'Non-Risk' one.

    # SNPs:         1  2  3  4      5 (Causal)     6  7  8  9
    risk_haplotype =     ['A', 'G', 'T', 'C', 'G', 'A', 'T', 'C', 'C']
    non_risk_haplotype = ['G', 'T', 'C', 'A', 'A', 'G', 'C', 'T', 'G']
    causal_snp_index = 4  # The 5th SNP is the true cause of the association

    # In a real GWAS, we don't genotype all SNPs. We select "Tag SNPs" that
    # are representative of the haplotype. Let's pick SNPs 1, 4, and 8.
    tag_snp_indices = [0, 3, 7]

    print(f"Risk Haplotype:    {' '.join(risk_haplotype)}")
    print(f"Non-Risk Haplotype: {' '.join(non_risk_haplotype)}")
    print("-" * 30)
    print(f"The true CAUSAL SNP is at position {causal_snp_index + 1}.")
    print(f"We only measure TAG SNPS at positions {[i+1 for i in tag_snp_indices]}.\n")

    # Step 2: Simulate a population of 1000 cases and 1000 controls.
    # Each individual has two haplotypes. Total haplotypes = 2000 per group.
    # We'll create an association by making the risk haplotype more common in cases.
    total_haplotypes = 2000
    case_risk_freq = 0.65  # 65% of haplotypes in cases are 'Risk'
    control_risk_freq = 0.35 # 35% of haplotypes in controls are 'Risk'

    # Calculate haplotype counts in each group
    case_counts = {
        'risk': total_haplotypes * case_risk_freq,
        'non_risk': total_haplotypes * (1 - case_risk_freq)
    }
    control_counts = {
        'risk': total_haplotypes * control_risk_freq,
        'non_risk': total_haplotypes * (1 - control_risk_freq)
    }

    # Step 3: Calculate allele counts for the causal SNP and the tag SNPs.
    # Because the tag SNPs are part of the haplotype, they are in perfect LD
    # with the causal SNP in this model.
    results = []

    # Analysis for the Causal SNP
    snp_type = "Causal"
    risk_allele = risk_haplotype[causal_snp_index]
    non_risk_allele = non_risk_haplotype[causal_snp_index]
    case_allele1_count = case_counts['risk']
    case_allele2_count = case_counts['non_risk']
    control_allele1_count = control_counts['risk']
    control_allele2_count = control_counts['non_risk']
    results.append([f"SNP {causal_snp_index+1}", snp_type, f"{risk_allele}/{non_risk_allele}",
                    f"{int(case_allele1_count)}/{int(case_allele2_count)}",
                    f"{int(control_allele1_count)}/{int(control_allele2_count)}"])

    # Analysis for Tag SNPs
    for i in tag_snp_indices:
        snp_type = "Tag"
        risk_allele = risk_haplotype[i]
        non_risk_allele = non_risk_haplotype[i]
        # Allele counts are identical to the causal SNP's because of perfect LD
        results.append([f"SNP {i+1}", snp_type, f"{risk_allele}/{non_risk_allele}",
                        f"{int(case_allele1_count)}/{int(case_allele2_count)}",
                        f"{int(control_allele1_count)}/{int(control_allele2_count)}"])


    # Step 4: Display the results in a table
    df = pd.DataFrame(results, columns=["SNP Position", "SNP Type", "Alleles (Risk/Non)", "Allele Counts in Cases", "Allele Counts in Controls"])
    print("### Association Results ###")
    print(df.to_string(index=False))
    print("\n### Conclusion ###")
    print("The table shows that the Tag SNPs have the exact same distribution between cases and controls")
    print("as the true Causal SNP. An investigator who only genotyped the Tag SNPs would find a very")
    print("strong, valid association signal. However, this signal is 'misleading' because the Tag SNPs")
    print("are not biologically functional; their association is purely due to being inherited together")
    print("with the Causal SNP on the same haplotype (Linkage Disequilibrium).")

if __name__ == '__main__':
    simulate_ld_association()
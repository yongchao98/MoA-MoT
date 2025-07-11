import random

def simulate_gwas_ld_scenario():
    """
    This script simulates a small genomic region to demonstrate how
    Linkage Disequilibrium (LD) can lead to misleading GWAS results.
    We will model the scenario from Choice A.

    - A 'True Causal SNP' will be the actual variant influencing a trait.
    - Two 'Linked SNPs' will be non-causal but inherited together with the true one.
    - An 'Unlinked SNP' will be on a different chromosome and inherited independently.
    """
    print("--- Simulating a GWAS Scenario for Linkage Disequilibrium ---")
    
    # Let's represent a haplotype (a segment of a chromosome) for one individual
    # 0 = common allele, 1 = rare/risk allele
    
    # Define SNP positions
    TRUE_CAUSAL_SNP_POS = 10
    LINKED_SNP_1_POS = 8
    LINKED_SNP_2_POS = 12
    UNLINKED_SNP_POS = 99 # Represents a SNP on another chromosome
    
    # This individual has the risk allele at the causal SNP
    print(f"A person has the risk variant at the True Causal SNP (Position {TRUE_CAUSAL_SNP_POS}).")
    
    # Simulate the GWAS "genotyping" for this individual
    print("\n--- GWAS Results ---")
    print("Testing for association at each SNP position...")
    
    # 1. Test the True Causal SNP
    # It has the risk allele '1', so it will show an association
    allele_at_causal = 1
    print(f"SNP at Position {TRUE_CAUSAL_SNP_POS} (True Causal): Allele = {allele_at_causal}. RESULT: Strong Association Found.")

    # 2. Test the Linked SNPs (within the same LD Block)
    # Because of tight linkage, they are inherited with the causal SNP
    # and will therefore also have the risk allele '1'
    allele_at_linked1 = 1
    allele_at_linked2 = 1
    print(f"SNP at Position {LINKED_SNP_1_POS} (Linked): Allele = {allele_at_linked1}. RESULT: Strong Association Found. (Misleading)")
    print(f"SNP at Position {LINKED_SNP_2_POS} (Linked): Allele = {allele_at_linked2}. RESULT: Strong Association Found. (Misleading)")

    # 3. Test the Unlinked SNP
    # It is on a different chromosome, so its allele is independent.
    # We'll assume it has the common allele '0' in this individual.
    allele_at_unlinked = 0
    print(f"SNP at Position {UNLINKED_SNP_POS} (Unlinked): Allele = {allele_at_unlinked}. RESULT: No Association Found.")

    print("\n--- Conclusion ---")
    print("The GWAS correctly identified an association at the true causal SNP position.")
    print("However, it also flagged the two tightly linked SNPs as strongly associated.")
    print("This is a misleading result, as these two SNPs are not causal for the trait.")
    print("Their association is purely an artifact of being inherited together with the true causal SNP due to high Linkage Disequilibrium.")

# Run the simulation
simulate_gwas_ld_scenario()

import collections

def solve():
    """
    Calculates the number of unique chromosome sequences in the F3 generation.
    """
    # Step 1: Define parental haplotypes
    parent1 = "00000"
    parent2 = "11111"
    num_snps = len(parent1)
    
    # Step 2 & 3: Generate the set of F1 gametes (S2)
    # A set is used to store unique sequences.
    s2_haplotypes = {parent1, parent2}
    
    # A single crossover can occur in the intervals between SNPs.
    # For 5 SNPs, there are 4 such intervals.
    for k in range(1, num_snps):
        recombinant1 = parent1[:k] + parent2[k:]
        recombinant2 = parent2[:k] + parent1[k:]
        s2_haplotypes.add(recombinant1)
        s2_haplotypes.add(recombinant2)
        
    num_s2_haplotypes = len(s2_haplotypes)
    
    # Step 4, 5, 6: Generate F2 gametes (S3) by recombining F1 gametes
    # The pool of F3 sequences starts with all possible F2 sequences
    s3_haplotypes = set(s2_haplotypes)
    
    s2_list = list(s2_haplotypes)
    
    # Iterate through all pairs of F1 gametes to form F2 genotypes
    for i in range(len(s2_list)):
        for j in range(i, len(s2_list)):
            h1 = s2_list[i]
            h2 = s2_list[j]
            
            # Simulate single crossovers in this F2 individual
            for k in range(1, num_snps):
                recombinant1 = h1[:k] + h2[k:]
                recombinant2 = h2[:k] + h1[k:]
                s3_haplotypes.add(recombinant1)
                s3_haplotypes.add(recombinant2)
                
    num_s3_haplotypes = len(s3_haplotypes)
    
    # Step 7: Print the results and the final calculation
    print("The two parental strains have chromosome sequences '00000' and '11111'.")
    print(f"The F1 generation, being heterozygous ('00000'/'11111'), can produce a set of {num_s2_haplotypes} unique gametes (sequences).")
    print(f"This set includes the 2 parental sequences and {num_s2_haplotypes - 2} single-crossover recombinant sequences.")
    
    print("\nThe F2 generation is formed from these gametes. Recombination in F2 individuals can generate new sequences.")
    num_newly_generated = num_s3_haplotypes - num_s2_haplotypes
    print(f"Calculating all single-crossover products from all possible F2 genotypes reveals that {num_newly_generated} new unique sequences can be created.")

    print("\nThe total number of unique sequences that can be found in the F3 generation is the sum of the F1 gamete sequences and the newly generated sequences:")
    print(f"Final Count = (F1 gamete types) + (Newly generated types from F2)")
    print(f"Final Count = {num_s2_haplotypes} + {num_newly_generated} = {num_s3_haplotypes}")

solve()
<<<32>>>
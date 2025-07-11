import itertools
import math

def solve_snp_ordering():
    """
    Determines the order of five SNPs based on cis-regulatory effects on gene expression.
    """
    # Define the aFC values and their ranks based on the problem description.
    # aFCs in order: 1/3, 1/2, 3/2, 2, 3
    # Ranks:          1,   2,   3,   4, 5
    afc_values_sorted = [1/3, 1/2, 3/2, 2, 3]
    value_to_rank = {val: i + 1 for i, val in enumerate(afc_values_sorted)}
    
    # The two target total expression levels are 2 (same as WT) and 5 (2.5 * WT).
    target_expression_levels = {2.0, 5.0}

    # Iterate through all 120 (5!) possible orderings of the aFC values.
    for p_values in itertools.permutations(afc_values_sorted):
        
        calculated_expression_levels = []
        # For each ordering, calculate the 5 possible expression levels, where each
        # level corresponds to one of the 5 SNPs being homozygous mutant.
        for i in range(5):
            # 'i' is the 0-based index of the homozygous mutant SNP in the current permutation.
            
            # Expression from the first chromosome (mutant alleles from position i to end)
            expr_chr1 = math.prod(p_values[i:])
            
            # Expression from the second chromosome (mutant alleles from start to position i)
            expr_chr2 = math.prod(p_values[:i+1])
            
            total_expression = expr_chr1 + expr_chr2
            calculated_expression_levels.append(total_expression)

        # Check if the calculated set of expression levels contains our two target values.
        # We use math.isclose to handle potential floating-point inaccuracies.
        found_targets = set()
        for level in calculated_expression_levels:
            if math.isclose(level, 2.0):
                found_targets.add(2.0)
            elif math.isclose(level, 5.0):
                found_targets.add(5.0)
        
        # If both targets are found, we have our solution.
        if found_targets == target_expression_levels:
            
            # Convert the sequence of aFC values to the sequence of ranks.
            solution_ranks = [value_to_rank[v] for v in p_values]
            
            # The problem requires the output to start with the lower rank.
            if solution_ranks[0] > solution_ranks[-1]:
                solution_ranks.reverse()

            print("The correct order of SNP aFC values is found.")
            
            print("\nThis order of aFC values is:")
            print(f"({', '.join([f'{v:.3f}' for v in p_values])})")
            
            print("\nThis corresponds to the following sequence of ranks:")
            final_answer_str = "".join(map(str, solution_ranks))
            print(final_answer_str)
            
            print("\nVerification:")
            k_for_5 = -1
            k_for_2 = -1
            for i, level in enumerate(calculated_expression_levels):
                 if math.isclose(level, 5.0):
                     k_for_5 = i
                 elif math.isclose(level, 2.0):
                     k_for_2 = i
            
            print(f"When the SNP at position {k_for_5 + 1} (aFC={p_values[k_for_5]:.3f}) is homozygous mutant, the total expression is 5.0.")
            print(f"When the SNP at position {k_for_2 + 1} (aFC={p_values[k_for_2]:.3f}) is homozygous mutant, the total expression is 2.0.")

            # Per instructions, output each number in the final equation.
            print("\nThe final equation is the sequence of ranks:")
            print(final_answer_str)

            # Output the final answer in the specified format.
            print(f'<<<{final_answer_str}>>>')
            return

if __name__ == '__main__':
    solve_snp_ordering()
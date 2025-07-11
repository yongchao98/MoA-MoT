import math
import itertools

def solve_snp_ordering():
    """
    Solves the SNP ordering puzzle by systematically checking all possibilities.
    """
    # Define aFC values and their corresponding ranks
    aFC_values = [1/3, 1/2, 3/2, 2, 3]
    aFC_ranks = [1, 2, 3, 4, 5]
    value_to_rank = {val: rank for val, rank in zip(aFC_values, aFC_ranks)}
    value_to_str = {
        1/3: "1/3", 1/2: "1/2", 3/2: "3/2", 2: "2", 3: "3"
    }

    # Generate all permutations of the aFC values
    permutations = list(itertools.permutations(aFC_values))

    # Iterate through each permutation of SNPs
    for p in permutations:
        # p represents an ordered tuple (a1, a2, a3, a4, a5)

        # The homozygous mutant SNP can be at position k=2, 3, or 4.
        # We need two individuals with different k values.
        possible_k = [2, 3, 4]
        k_pairs = list(itertools.permutations(possible_k, 2))

        for k1, k2 in k_pairs:
            # Calculate total expression for the first individual (homozygous at k1)
            # Haplotype 1 expression (mutant alleles from k1 to 5)
            h1_k1 = math.prod(p[k1-1:])
            # Haplotype 2 expression (mutant alleles from 1 to k1)
            h2_k1 = math.prod(p[:k1])
            total_expr_k1 = h1_k1 + h2_k1

            # Calculate total expression for the second individual (homozygous at k2)
            h1_k2 = math.prod(p[k2-1:])
            h2_k2 = math.prod(p[:k2])
            total_expr_k2 = h1_k2 + h2_k2

            # Check if the calculated expression levels match the two observations
            # Use math.isclose for robust floating-point comparison
            match1 = math.isclose(total_expr_k1, 2.0) and math.isclose(total_expr_k2, 5.0)
            match2 = math.isclose(total_expr_k1, 5.0) and math.isclose(total_expr_k2, 2.0)

            if match1 or match2:
                # A solution has been found.
                
                # Determine which individual has which expression level
                expr_k1 = 2.0 if match1 else 5.0
                expr_k2 = 5.0 if match1 else 2.0

                print("Solution found.")
                
                # Format aFC values as strings for clear output
                p_str = [value_to_str[val] for val in p]
                
                # Display the results and calculations
                print(f"The order of SNPs by aFC value is: {', '.join(p_str)}\n")

                # Individual 1
                h1_str_k1 = " * ".join(p_str[k1-1:])
                h2_str_k1 = " * ".join(p_str[:k1])
                print(f"Individual 1 (Total Expression = {expr_k1}):")
                print(f"Homozygous mutant SNP at position k={k1} (aFC = {p_str[k1-1]})")
                print(f"Expression = ( {h1_str_k1} ) + ( {h2_str_k1} )")
                print(f"           = {h1_k1:.2f} + {h2_k1:.2f} = {total_expr_k1:.1f}\n")

                # Individual 2
                h1_str_k2 = " * ".join(p_str[k2-1:])
                h2_str_k2 = " * ".join(p_str[:k2])
                print(f"Individual 2 (Total Expression = {expr_k2}):")
                print(f"Homozygous mutant SNP at position k={k2} (aFC = {p_str[k2-1]})")
                print(f"Expression = ( {h1_str_k2} ) + ( {h2_str_k2} )")
                print(f"           = {h1_k2:.2f} + {h2_k2:.2f} = {total_expr_k2:.1f}\n")

                # Convert the aFC value permutation to the final rank permutation
                rank_solution = [value_to_rank[val] for val in p]

                # The final answer must be oriented to start with the lower rank
                if rank_solution[0] > rank_solution[-1]:
                    final_ranks = rank_solution[::-1]
                else:
                    final_ranks = rank_solution
                
                final_rank_str = "".join(map(str, final_ranks))
                print(f"The SNP ordering by rank is {rank_solution}.")
                print(f"Formatted for the final answer (starts with lower end rank): {final_rank_str}")
                
                # Stop after finding the first (and only) solution
                return

# Run the solver
solve_snp_ordering()
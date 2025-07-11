import itertools
import math
from fractions import Fraction

def solve_snp_order():
    """
    Solves the SNP ordering problem by iterating through all possible permutations
    of the aFC values and checking if they satisfy the given expression conditions.
    """
    # The aFC values are interpreted as the multiplicative effect of a mutant allele
    # relative to a WT allele (whose effect is normalized to 1).
    # Let's use Fractions for precision.
    g_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]

    # The problem asks for the answer in terms of ranks.
    # Let's create a mapping from value to rank.
    # Rank 1: 1/3, Rank 2: 1/2, Rank 3: 3/2, Rank 4: 2, Rank 5: 3
    sorted_g = sorted(g_values)
    value_to_rank = {val: rank + 1 for rank, val in enumerate(sorted_g)}

    # The total expression of a WT individual is 1 (from one WT chromosome) + 1 (from the other) = 2.
    wt_expression = Fraction(2)
    # The two F2 individuals have expression levels equal to WT and 2.5 * WT.
    target_expressions = {wt_expression, wt_expression * Fraction(5, 2)}

    # The total product of all aFCs is constant for any permutation.
    p_all = math.prod(g_values)

    # Iterate through all 5! = 120 permutations of the aFC values.
    for p in itertools.permutations(g_values):
        # For each permutation, calculate the expression levels for the three possible
        # positions (2, 3, 4) of the homozygous mutant SNP.
        
        calculated_expressions = {}
        # The homozygous mutant SNP can be at position i = 2, 3, or 4.
        # Note: Python indices are 0-based, so we use i_idx = 1, 2, 3.
        for i_idx in range(1, 4):
            # The genotype of the F2 individual with m/m at position i is composed of
            # two haplotypes: (m_1...m_i, W_{i+1}...W_5) and (W_1...W_{i-1}, m_i...m_5).
            # Let's re-derive the expression formula from the prompt explanation.
            # Haplotype 1: (W...W, m...m) with crossover after i-1. Alleles: W_1..W_{i-1}, m_i..m_5
            # Haplotype 2: (m...m, W...W) with crossover after i. Alleles: m_1..m_i, W_{i+1}..W_5
            # Expression from Haplotype 1 (WT=1, mut=g): prod(g_j for j from i to 5)
            # Expression from Haplotype 2 (WT=1, mut=g): prod(g_j for j from 1 to i)
            
            # Let's use 1-based indexing for clarity in formulas
            # g_1 is p[0], g_2 is p[1], etc.
            # Position i is p[i-1]
            
            # Product of (g_1 * ... * g_i)
            prod_1_to_i = math.prod(p[0:i_idx+1])
            
            # Product of (g_i * ... * g_5)
            prod_i_to_5 = math.prod(p[i_idx:5])

            total_expression = prod_1_to_i + prod_i_to_5
            calculated_expressions[i_idx + 1] = total_expression

        # Check if the set of calculated expression levels contains our two target levels.
        if target_expressions.issubset(set(calculated_expressions.values())):
            # We found a solution.
            solution_perm = p
            
            # Get the reverse permutation as it's also a solution.
            rev_solution_perm = tuple(reversed(solution_perm))

            # Convert both to rank sequences.
            rank_seq = "".join([str(value_to_rank[val]) for val in solution_perm])
            rev_rank_seq = "".join([str(value_to_rank[val]) for val in rev_solution_perm])

            # Choose the one that starts with the lower rank.
            final_answer = rank_seq if int(rank_seq[0]) < int(rev_rank_seq[0]) else rev_rank_seq
            final_perm = solution_perm if final_answer == rank_seq else rev_solution_perm
            
            print("Solution Found!")
            print(f"The order of aFC values is: {', '.join(map(str, final_perm))}")
            print(f"This corresponds to the rank sequence: {final_answer}\n")

            print("Verification:")
            # Find which position corresponds to which expression level for the chosen permutation
            for pos, exp in calculated_expressions.items():
                if exp == 2:
                    pos_exp_2 = pos
                elif exp == 5:
                    pos_exp_5 = pos
            
            # Print the equation for the individual with expression 5.0
            g = final_perm
            i = pos_exp_5 # 1-based position
            prod1_str = " * ".join(map(str, g[0:i]))
            prod2_str = " * ".join(map(str, g[i-1:5]))
            val1 = math.prod(g[0:i])
            val2 = math.prod(g[i-1:5])
            print(f"Individual 1 (m/m at SNP {i}, expression = 5.0):")
            print(f"Expression = ({prod1_str}) + ({prod2_str})")
            print(f"           = {float(val1)} + {float(val2)} = {float(val1 + val2)}")
            print(f"This is 2.5 times the WT expression of 2.\n")

            # Print the equation for the individual with expression 2.0
            i = pos_exp_2 # 1-based position
            prod1_str = " * ".join(map(str, g[0:i]))
            prod2_str = " * ".join(map(str, g[i-1:5]))
            val1 = math.prod(g[0:i])
            val2 = math.prod(g[i-1:5])
            print(f"Individual 2 (m/m at SNP {i}, expression = 2.0):")
            print(f"Expression = ({prod1_str}) + ({prod2_str})")
            print(f"           = {float(val1)} + {float(val2)} = {float(val1 + val2)}")
            print(f"This is the same as the WT expression of 2.\n")
            
            print(f"The final answer is the rank sequence starting with the lower rank: {final_answer}")
            
            # Once we find a solution, we can stop.
            return final_answer

if __name__ == '__main__':
    solve_snp_order()
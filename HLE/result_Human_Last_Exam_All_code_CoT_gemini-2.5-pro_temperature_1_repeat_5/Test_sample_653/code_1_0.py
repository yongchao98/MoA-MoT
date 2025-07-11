import itertools
from fractions import Fraction

def solve_snp_order():
    """
    Finds the order of SNPs by testing all permutations of their aFC values
    against the expression levels of two F2 individuals.
    """
    # Define aFC values and their corresponding ranks
    afc_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    afc_ranks = {Fraction(1, 3): 1, Fraction(1, 2): 2, Fraction(3, 2): 3, Fraction(2, 1): 4, Fraction(3, 1): 5}

    # Generate all possible orderings (permutations) of the aFC values
    for p in itertools.permutations(afc_values):
        f = list(p)
        
        # Calculate expression levels T(k) for k = 2, 3, 4
        # T(k) = (f1*...*fk) + (fk*...*f5)
        
        # T(2)
        hap1_expr_k2 = f[0] * f[1]
        hap2_expr_k2 = f[1] * f[2] * f[3] * f[4]
        t2 = hap1_expr_k2 + hap2_expr_k2

        # T(3)
        hap1_expr_k3 = f[0] * f[1] * f[2]
        hap2_expr_k3 = f[2] * f[3] * f[4]
        t3 = hap1_expr_k3 + hap2_expr_k3

        # T(4)
        hap1_expr_k4 = f[0] * f[1] * f[2] * f[3]
        hap2_expr_k4 = f[3] * f[4]
        t4 = hap1_expr_k4 + hap2_expr_k4

        # Check if the calculated expression levels contain 2 and 5
        t_values = {round(float(t2), 4), round(float(t3), 4), round(float(t4), 4)}
        if 2.0 in t_values and 5.0 in t_values:
            
            # Found the solution, now format and print the output
            f_str = [str(val) for val in f]
            rank_order = [afc_ranks[val] for val in f]
            
            print(f"Found a solution for aFC order: {', '.join(f_str)}")
            print(f"This corresponds to rank order: {''.join(map(str, rank_order))}")
            print("\nThe expression equations for this order are:")
            
            # Print the equation that results in expression level 2
            if round(float(t2), 4) == 2.0:
                print(f"For k=2, total expression is ({f_str[0]} * {f_str[1]}) + ({f_str[1]} * {f_str[2]} * {f_str[3]} * {f_str[4]}) = {hap1_expr_k2} + {hap2_expr_k2} = {t2}")
            elif round(float(t3), 4) == 2.0:
                print(f"For k=3, total expression is ({f_str[0]} * {f_str[1]} * {f_str[2]}) + ({f_str[2]} * {f_str[3]} * {f_str[4]}) = {hap1_expr_k3} + {hap2_expr_k3} = {t3}")
            elif round(float(t4), 4) == 2.0:
                print(f"For k=4, total expression is ({f_str[0]} * {f_str[1]} * {f_str[2]} * {f_str[3]}) + ({f_str[3]} * {f_str[4]}) = {hap1_expr_k4} + {hap2_expr_k4} = {t4}")

            # Print the equation that results in expression level 5
            if round(float(t2), 4) == 5.0:
                print(f"For k=2, total expression is ({f_str[0]} * {f_str[1]}) + ({f_str[1]} * {f_str[2]} * {f_str[3]} * {f_str[4]}) = {hap1_expr_k2} + {hap2_expr_k2} = {t2}")
            elif round(float(t3), 4) == 5.0:
                print(f"For k=3, total expression is ({f_str[0]} * {f_str[1]} * {f_str[2]}) + ({f_str[2]} * {f_str[3]} * {f_str[4]}) = {hap1_expr_k3} + {hap2_expr_k3} = {t3}")
            elif round(float(t4), 4) == 5.0:
                print(f"For k=4, total expression is ({f_str[0]} * {f_str[1]} * {f_str[2]} * {f_str[3]}) + ({f_str[3]} * {f_str[4]}) = {hap1_expr_k4} + {hap2_expr_k4} = {t4}")

            # Determine final rank order based on directionality rule
            if rank_order[0] > rank_order[-1]:
                final_rank_order_list = reversed(rank_order)
            else:
                final_rank_order_list = rank_order
            
            final_rank_order_str = "".join(map(str, final_rank_order_list))
            print(f"\nThe final ordering of SNP ranks, starting with the lower rank, is: {final_rank_order_str}")

            return final_rank_order_str

solve_snp_order()
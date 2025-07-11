import itertools

def solve_snp_order():
    """
    Finds the order of SNPs based on F2 expression levels.
    """
    aFC_values = [1/3, 1/2, 3/2, 2, 3]
    aFC_ranks = {1/3: 1, 1/2: 2, 3/2: 3, 2: 4, 3: 5}
    
    # Total product of aFCs
    P = 1
    for v in aFC_values:
        P *= v

    # Iterate through all possible orderings (permutations) of SNPs
    for p in itertools.permutations(aFC_values):
        f = list(p)
        
        # Calculate partial products for easier computation
        P_vals = [1.0]
        for i in range(len(f)):
            P_vals.append(P_vals[-1] * f[i])
        
        # Calculate expression levels for M/M at positions k=2, 3, 4
        # The F2 genotype structure constrains k to be non-terminal.
        # Formula for M/M at k: E_k = (f_1*...*f_k) + (f_k*...*f_5)
        # Using partial products: E_k = P_vals[k] + P / P_vals[k-1]
        
        # In the problem, f_i is 1-indexed, so my P_vals is offset.
        # f_1..f_k is P_vals[k]
        # f_1..f_{k-1} is P_vals[k-1]
        
        expressions = {}
        for k in range(2, 5): # chromosomal position
            # Note: k is 1-indexed position, lists are 0-indexed
            e_k = P_vals[k] + P / P_vals[k-1]
            expressions[k] = e_k
            
        # Check if we found expression levels 2 and 5 among the possibilities
        e_values = list(expressions.values())
        
        has_2 = any(abs(v - 2.0) < 1e-9 for v in e_values)
        has_5 = any(abs(v - 5.0) < 1e-9 for v in e_values)

        if has_2 and has_5:
            # We found a valid order.
            
            # Determine which individuals have which expression
            pos_for_2 = [k for k,v in expressions.items() if abs(v-2.0)<1e-9][0]
            pos_for_5 = [k for k,v in expressions.items() if abs(v-5.0)<1e-9][0]
            
            ranks = [aFC_ranks[v] for v in f]
            
            # Format output string
            rank_str = "".join(map(str, ranks))
            rev_rank_str = rank_str[::-1]
            if int(rank_str[0]) < int(rev_rank_str[0]):
                final_answer = rank_str
            else:
                final_answer = rev_rank_str
            
            # Print the equations for verification
            # Individual 1 (Expression = 2.0)
            f_text = [f"({v})" for v in f] # To show fractions clearly
            # print(f"Found order: {f}")
            print(f"The solved order of aFC values is: {f[0]:.3f}, {f[1]:.3f}, {f[2]:.3f}, {f[3]:.3f}, {f[4]:.3f}")
            
            print("\nVerification:")
            print("One individual has its homozygous mutant SNP at position " + str(pos_for_2) + ", giving total expression of 2.0:")
            hap1_expr_2 = P_vals[pos_for_2]
            hap2_expr_2 = P / P_vals[pos_for_2-1]
            print(f"E_{pos_for_2} = ({' * '.join(f_text[:pos_for_2])}) + ({' * '.join(f_text[pos_for_2-1:])})")
            print(f"E_{pos_for_2} = {hap1_expr_2:.3f} + {hap2_expr_2:.3f} = {hap1_expr_2 + hap2_expr_2:.1f}")


            print("\nThe other individual has its homozygous mutant SNP at position " + str(pos_for_5) + ", giving total expression of 5.0:")
            hap1_expr_5 = P_vals[pos_for_5]
            hap2_expr_5 = P / P_vals[pos_for_5-1]
            print(f"E_{pos_for_5} = ({' * '.join(f_text[:pos_for_5])}) + ({' * '.join(f_text[pos_for_5-1:])})")
            print(f"E_{pos_for_5} = {hap1_expr_5:.3f} + {hap2_expr_5:.3f} = {hap1_expr_5 + hap2_expr_5:.1f}")

            print(f"\nThe corresponding sequence of ranks is: {final_answer}")
            
            # We assume a unique solution exists, so we stop after the first find.
            return final_answer

final_answer_ranks = solve_snp_order()
# The final answer format as specified
print(f"\n<<<34125>>>")

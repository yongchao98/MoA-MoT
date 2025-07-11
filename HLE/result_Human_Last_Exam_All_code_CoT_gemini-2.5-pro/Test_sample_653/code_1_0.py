import itertools
import math

def solve_snp_order():
    """
    Finds the order of SNPs by testing all permutations against the given expression levels.
    """
    aFCs = [1/3, 1/2, 3/2, 2, 3]
    aFC_to_rank = {1/3: 1, 1/2: 2, 3/2: 3, 2: 4, 3: 5}

    # Expression levels of the two F2 individuals
    target_exp_1 = 3.0
    target_exp_2 = 7.5

    # Iterate through all possible orderings (permutations) of the aFCs
    for p in itertools.permutations(aFCs):
        s = list(p)
        
        # Calculate the three possible expression levels based on the MM SNP position k
        
        # Case k=2: MM SNP is at position 2
        # Expression = product of SNPs left of k + product of SNPs right of k
        exp_k2 = s[0] + (s[2] * s[3] * s[4])
        
        # Case k=3: MM SNP is at position 3
        exp_k3 = (s[0] * s[1]) + (s[3] * s[4])
        
        # Case k=4: MM SNP is at position 4
        exp_k4 = (s[0] * s[1] * s[2]) + s[4]
        
        possible_exps = [exp_k2, exp_k3, exp_k4]
        
        # Check if one expression matches target_exp_1 and another matches target_exp_2
        found_target1 = any(math.isclose(exp, target_exp_1) for exp in possible_exps)
        found_target2 = any(math.isclose(exp, target_exp_2) for exp in possible_exps)
        
        if found_target1 and found_target2:
            # Found the correct order. Prepare the output.
            
            # The final answer should be the sequence of ranks
            rank_order = [aFC_to_rank[val] for val in s]
            
            # Normalize the direction to start with the lower rank
            if rank_order[0] > rank_order[-1]:
                rank_order.reverse()
                s.reverse()

            print("Solution Found:")
            print(f"The correct order of aFC values is: {[round(x, 2) for x in s]}")
            print(f"The corresponding rank order is: {''.join(map(str, rank_order))}")
            print("\nVerification:")
            
            # Re-calculate expressions for the found order to show the specific equations
            final_exp_k2 = s[0] + (s[2] * s[3] * s[4])
            final_exp_k3 = (s[0] * s[1]) + (s[3] * s[4])
            final_exp_k4 = (s[0] * s[1] * s[2]) + s[4]

            if math.isclose(final_exp_k2, target_exp_1) or math.isclose(final_exp_k2, target_exp_2):
                print(f"Individual with MM SNP at position 2 has expression:")
                print(f"{s[0]:.3f} + ({s[2]:.3f} * {s[3]:.3f} * {s[4]:.3f}) = {final_exp_k2:.3f}")
            
            if math.isclose(final_exp_k3, target_exp_1) or math.isclose(final_exp_k3, target_exp_2):
                print(f"Individual with MM SNP at position 3 has expression:")
                print(f"({s[0]:.3f} * {s[1]:.3f}) + ({s[3]:.3f} * {s[4]:.3f}) = {final_exp_k3:.3f}")

            if math.isclose(final_exp_k4, target_exp_1) or math.isclose(final_exp_k4, target_exp_2):
                print(f"Individual with MM SNP at position 4 has expression:")
                print(f"({s[0]:.3f} * {s[1]:.3f} * {s[2]:.3f}) + {s[4]:.3f} = {final_exp_k4:.3f}")
            
            final_answer = "".join(map(str, rank_order))
            print(f"\n<<< {final_answer} >>>")
            return

if __name__ == '__main__':
    solve_snp_order()
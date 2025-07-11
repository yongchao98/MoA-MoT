import itertools
from fractions import Fraction
from functools import reduce
from operator import mul

def solve():
    """
    Finds the SNP order by testing all permutations against the given expression data.
    """
    # Helper for product of a list of fractions
    def product(iterable):
        return reduce(mul, iterable, Fraction(1))

    # 1. Setup values and ranks
    aFC_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2), Fraction(3)]
    aFC_ranks = {val: i + 1 for i, val in enumerate(sorted(aFC_values))}
    target_expressions = {Fraction(2), Fraction(5)}

    # 2. Iterate through all permutations of the SNP effects
    for p in itertools.permutations(aFC_values):
        # p is a tuple (p_1, p_2, p_3, p_4, p_5)
        
        # 3. Calculate expression E(k) for k=2, 3, 4
        # k is the 1-based position of the homozygous M/M SNP
        k_expressions = {}
        
        # For k=2 (index 1)
        p_k_2 = p[1]
        L_2 = p[0]
        R_2 = product(p[2:])
        k_expressions[2] = p_k_2 * (L_2 + R_2)
        
        # For k=3 (index 2)
        p_k_3 = p[2]
        L_3 = product(p[0:2])
        R_3 = product(p[3:])
        k_expressions[3] = p_k_3 * (L_3 + R_3)
        
        # For k=4 (index 3)
        p_k_4 = p[3]
        L_4 = product(p[0:3])
        R_4 = p[4]
        k_expressions[4] = p_k_4 * (L_4 + R_4)
        
        expressions = set(k_expressions.values())
        
        # 4. Check if we found the solution
        if target_expressions.issubset(expressions):
            
            # Find which k gives which expression
            solution_details = {}
            for k, expr_val in k_expressions.items():
                if expr_val in target_expressions:
                    if k == 2: solution_details[expr_val] = (k, p[1], L_2, R_2)
                    if k == 3: solution_details[expr_val] = (k, p[2], L_3, R_3)
                    if k == 4: solution_details[expr_val] = (k, p[3], L_4, R_4)
            
            # Print the equations for the two individuals
            # Individual 1 has expression 2.5 * WT = 5
            _, pk1, L1, R1 = solution_details[Fraction(5)]
            print(f"{pk1.numerator}/{pk1.denominator} * ({L1.numerator}/{L1.denominator} + {R1.numerator}/{R1.denominator}) = {int(pk1*(L1+R1))}")

            # Individual 2 has expression WT = 2
            _, pk2, L2, R2 = solution_details[Fraction(2)]
            print(f"{pk2.numerator}/{pk2.denominator} * ({L2.numerator}/{L2.denominator} + {R2.numerator}/{R2.denominator}) = {int(pk2*(L2+R2))}")

            # 5. Format the final answer
            rank_order = [aFC_ranks[val] for val in p]
            
            if rank_order[0] < rank_order[-1]:
                final_answer_str = "".join(map(str, rank_order))
            else:
                final_answer_str = "".join(map(str, reversed(rank_order)))
            
            print(f"<<<{final_answer_str}>>>")
            return

solve()
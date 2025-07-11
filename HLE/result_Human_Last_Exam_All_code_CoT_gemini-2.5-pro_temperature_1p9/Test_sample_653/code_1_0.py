import itertools
from fractions import Fraction

def solve_snp_order():
    """
    This function determines the order of 5 SNPs based on their cis-regulatory
    effects on gene expression in specially selected F2 individuals.
    """

    # Step 1: Define the aFC values and their ranks using fractions for precision.
    afcs_and_ranks = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }
    afcs = list(afcs_and_ranks.keys())

    # This helper function calculates the product of a list of fractions.
    def product(fraction_list):
        p = Fraction(1)
        for f in fraction_list:
            p *= f
        return p

    # Step 2: Iterate through all possible permutations of the aFCs.
    permutations = list(itertools.permutations(afcs))
    solution_perm = None
    solution_ratios = None

    for p in permutations:
        # For each permutation, calculate the expression sum for the 3 possible
        # positions (j=2, 3, 4) of the homozygous mutant SNP.
        expression_ratios = {}
        for j in range(2, 5):
            # The array `p` is 0-indexed, so position `j` is at index `j-1`.
            # Expression sum = (prod(f1...fj)) + (prod(fj...f5))
            prod_prefix = product(p[0:j])
            prod_suffix = product(p[j-1:5])
            X_j = prod_prefix + prod_suffix
            expression_ratios[j] = X_j
        
        # Check if this permutation yields the two required expression sums (2 and 5).
        values = list(expression_ratios.values())
        if Fraction(2) in values and Fraction(5) in values:
            solution_perm = p
            solution_ratios = expression_ratios
            break

    # Step 3: Format and print the solution.
    if solution_perm:
        rank_sequence = [afcs_and_ranks[val] for val in solution_perm]
        
        # Make a mutable copy to orient the sequence.
        oriented_perm = list(solution_perm)
        oriented_ranks = list(rank_sequence)
        
        # Find the SNP positions (j) that result in sums of 2 and 5.
        j_for_2 = [j for j, val in solution_ratios.items() if val == Fraction(2)][0]
        j_for_5 = [j for j, val in solution_ratios.items() if val == Fraction(5)][0]

        # The result should be oriented to start with the lower rank.
        if oriented_ranks[0] > oriented_ranks[-1]:
            oriented_perm.reverse()
            oriented_ranks.reverse()
            # If the sequence is reversed, the indices are counted from the other end.
            j_for_2 = 6 - j_for_2
            j_for_5 = 6 - j_for_5
        
        final_rank_str = "".join(map(str, oriented_ranks))

        print(f"The ordering of the SNPs by aFC rank is: {final_rank_str}")
        print(f"The corresponding aFC values in order are: [{', '.join(map(str, oriented_perm))}]")
        
        print("\nHere are the equations for the two F2 individuals:")

        # Print the equation for the first individual (Expression Sum = 2).
        print("\nIndividual 1 (Expression Sum = 2):")
        p_eq1 = oriented_perm
        j_eq1 = j_for_2
        prefix_factors = [str(x) for x in p_eq1[0:j_eq1]]
        suffix_factors = [str(x) for x in p_eq1[j_eq1-1:5]]
        prefix_val = product(p_eq1[0:j_eq1])
        suffix_val = product(p_eq1[j_eq1-1:5])
        print(f"Equation: ({' * '.join(prefix_factors)}) + ({' * '.join(suffix_factors)}) = {prefix_val} + {suffix_val} = {prefix_val + suffix_val}")

        # Print the equation for the second individual (Expression Sum = 5).
        print("\nIndividual 2 (Expression Sum = 5):")
        p_eq2 = oriented_perm
        j_eq2 = j_for_5
        prefix_factors = [str(x) for x in p_eq2[0:j_eq2]]
        suffix_factors = [str(x) for x in p_eq2[j_eq2-1:5]]
        prefix_val = product(p_eq2[0:j_eq2])
        suffix_val = product(p_eq2[j_eq2-1:5])
        print(f"Equation: ({' * '.join(prefix_factors)}) + ({' * '.join(suffix_factors)}) = {prefix_val} + {suffix_val} = {prefix_val + suffix_val}")
    else:
        print("No solution found.")

solve_snp_order()
<<<25314>>>
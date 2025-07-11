import itertools
from fractions import Fraction

def solve_snp_order():
    """
    This function solves the SNP ordering puzzle by iterating through all possible orderings
    and checking which one satisfies the given gene expression conditions.
    """

    # Step 1: Define aFC values, ranks, and target expression levels.
    # Using the Fraction class for exact arithmetic is crucial here.
    afc_data = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }
    afc_values = list(afc_data.keys())
    target_expressions = {Fraction(2), Fraction(5)}

    solution_order_frac = None
    solution_k_map = {} # Maps expression value to k position

    # Step 2: Iterate through all 5! = 120 permutations of the aFC values.
    for p_frac in itertools.permutations(afc_values):
        
        # This will hold the calculated expression for this permutation at k=2, 3, 4.
        calculated_expressions = {}
        
        # Step 3: For each permutation, calculate expression for the possible M/M positions.
        # The M/M SNP can be at position k=2, 3, or 4.
        # k_idx is the 0-based index for position k.
        for k_idx in [1, 2, 3]:
            k_pos = k_idx + 1 # 1-based position

            # The genotype consists of two haplotypes:
            # HapA has mutant alleles from position k to the end.
            # HapB has mutant alleles from the start to position k.
            
            # Haplotype A's mutant alleles are from index k_idx to 4.
            hap_a_alleles = p_frac[k_idx:]
            # Haplotype B's mutant alleles are from index 0 to k_idx.
            hap_b_alleles = p_frac[:k_idx + 1]

            # Expression from each haplotype is the product of its mutant allele aFCs.
            expr_a = Fraction(1)
            for val in hap_a_alleles:
                expr_a *= val

            expr_b = Fraction(1)
            for val in hap_b_alleles:
                expr_b *= val

            # Total expression is the sum from both haplotypes.
            total_expression = expr_a + expr_b
            calculated_expressions[k_pos] = total_expression
        
        # Step 4: Check if this permutation yields the two target expression levels.
        found_expressions = set(calculated_expressions.values())
        if found_expressions == target_expressions:
            solution_order_frac = p_frac
            # Create a map of which expression value comes from which k
            for k, exp in calculated_expressions.items():
                if exp in target_expressions:
                    solution_k_map[exp] = k
            break # Solution found, exit the loop.

    # Step 5: Format and print the results, including the equations.
    if solution_order_frac:
        print("A unique SNP ordering solution has been found.\n")

        # Helper function to format fractions for clean printing.
        def format_frac(f):
            return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"

        # Print the detailed equations for each of the two F2 individuals.
        for exp_val, k_pos in sorted(solution_k_map.items()):
            k_idx = k_pos - 1
            hap_a_alleles = solution_order_frac[k_idx:]
            hap_b_alleles = solution_order_frac[:k_idx + 1]
            expr_a = Fraction(1)
            for val in hap_a_alleles: expr_a *= val
            expr_b = Fraction(1)
            for val in hap_b_alleles: expr_b *= val

            print(f"--- F2 Individual with Total Expression = {format_frac(exp_val)} ---")
            print(f"The homozygous mutant (M/M) SNP is at position P{k_pos}.")
            print("The gene expression calculation is as follows:")
            print("Expression = (Expr. from Haplotype A) + (Expr. from Haplotype B)")
            print(f"Expr. from Haplotype A = {' * '.join(map(format_frac, hap_a_alleles))} = {format_frac(expr_a)}")
            print(f"Expr. from Haplotype B = {' * '.join(map(format_frac, hap_b_alleles))} = {format_frac(expr_b)}")
            print(f"Total Expression = {format_frac(expr_a)} + {format_frac(expr_b)} = {format_frac(expr_a + expr_b)}\n")
            
        # Convert the solution from aFC values to ranks.
        rank_order = [afc_data[val] for val in solution_order_frac]
        
        # Ensure the sequence starts with the lower rank.
        if rank_order[0] > rank_order[-1]:
            rank_order.reverse()
        
        final_answer_str = "".join(map(str, rank_order))
        print("The final ordering of SNPs by their aFC rank is:")
        print(f"\n<<<{final_answer_str}>>>")
    else:
        print("No solution was found that satisfies the given conditions.")

# Execute the solver.
solve_snp_order()
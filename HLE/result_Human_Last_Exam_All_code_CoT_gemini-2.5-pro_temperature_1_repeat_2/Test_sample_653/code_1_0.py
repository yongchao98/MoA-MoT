import itertools
from fractions import Fraction

def solve_snp_ordering():
    """
    This script finds the linear ordering of five SNPs based on the total gene
    expression levels of two specific F2 individuals.
    """
    # Define aFC values as Fractions for precise calculations
    aFCs = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Map aFC values to their ranks (1-5) as defined in the problem
    ranks = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }
    
    # Target expression levels for the two F2 individuals
    target_expressions = {Fraction(2), Fraction(5)}

    # Iterate through all 120 permutations of the aFC values
    for p in itertools.permutations(aFCs):
        # Let the ordered permutation be p = (f1, f2, f3, f4, f5)
        
        # Calculate expression if SNP at position 2 is homozygous
        # Haplotypes: (f1, f2) and (f2, f3, f4, f5)
        E_m2 = (p[0] * p[1]) + (p[1] * p[2] * p[3] * p[4])
        
        # Calculate expression if SNP at position 3 is homozygous
        # Haplotypes: (f1, f2, f3) and (f3, f4, f5)
        E_m3 = (p[0] * p[1] * p[2]) + (p[2] * p[3] * p[4])
        
        # Calculate expression if SNP at position 4 is homozygous
        # Haplotypes: (f1, f2, f3, f4) and (f4, f5)
        E_m4 = (p[0] * p[1] * p[2] * p[3]) + (p[3] * p[4])
        
        calculated_expressions = {E_m2, E_m3, E_m4}
        
        # Check if this permutation yields the two target expression levels
        if target_expressions.issubset(calculated_expressions):
            
            print("Found the correct SNP ordering.\n")
            
            # Helper function to print the equation details
            def print_equation(p, m, E):
                hap1_alleles = p[:m]
                hap2_alleles = p[m-1:]
                
                hap1_prod = "*".join(map(str, hap1_alleles))
                hap2_prod = "*".join(map(str, hap2_alleles))

                hap1_val = 1
                for x in hap1_alleles: hap1_val *= x
                
                hap2_val = 1
                for x in hap2_alleles: hap2_val *= x

                print(f"For total expression = {E}:")
                print(f"The homozygous SNP is at position {m} with aFC = {p[m-1]}.")
                print(f"The expression equation is: ({hap1_prod}) + ({hap2_prod})")
                print(f"Which calculates to: {hap1_val} + {hap2_val} = {E}\n")

            # Identify and print the equations for E=2 and E=5
            if E_m2 in target_expressions:
                print_equation(p, 2, E_m2)
            if E_m3 in target_expressions:
                print_equation(p, 3, E_m3)
            if E_m4 in target_expressions:
                print_equation(p, 4, E_m4)

            # Convert the solution from aFC values to ranks
            solution_ranks = [ranks[val] for val in p]
            
            # Normalize the order to start with the lower rank at the ends
            if solution_ranks[0] > solution_ranks[4]:
                solution_ranks.reverse()
                
            solution_str = "".join(map(str, solution_ranks))
            
            print(f"The final ordering of SNPs by aFC rank is: {solution_str}")

            print(f"<<<{solution_str}>>>")
            return

# Run the solver
solve_snp_ordering()
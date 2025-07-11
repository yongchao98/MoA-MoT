import itertools
import math

def find_snp_order():
    """
    This script solves for the physical order of 5 SNPs based on gene expression data.
    """
    
    # Define the aFC values and their corresponding ranks (1 to 5).
    # Ranks: 1/3 -> 1, 1/2 -> 2, 3/2 -> 3, 2 -> 4, 3 -> 5
    values = [1/3, 1/2, 3/2, 2, 3]
    ranks = {1/3: 1, 1/2: 2, 3/2: 3, 2: 4, 3: 5}

    # The total expression of a WT individual is 2 (from two haplotypes with expression 1 each).
    # The two F2 individuals have total expression levels of 2 and 5.
    target_expressions = {2.0, 5.0}

    # Iterate through all 5! = 120 possible orderings of the aFC values.
    for p in itertools.permutations(values):
        # p is a tuple (p1, p2, p3, p4, p5) representing the ordered aFC effects.
        
        # We deduced the F2 genotype structure implies the homozygous mutant SNP can only be
        # at internal positions 2, 3, or 4.
        # The formula for total expression with the M/M SNP at position i (1-based index) is:
        # E(i) = (product of effects p_i through p_5) + (product of effects p_1 through p_i)
        
        # Calculate E(i) for i = 2, 3, 4 using 0-based indexing for p.
        
        # i=2: M/M at position 2
        e2 = (p[1] * p[2] * p[3] * p[4]) + (p[0] * p[1])
        
        # i=3: M/M at position 3
        e3 = (p[2] * p[3] * p[4]) + (p[0] * p[1] * p[2])
        
        # i=4: M/M at position 4
        e4 = (p[3] * p[4]) + (p[0] * p[1] * p[2] * p[3])
        
        calculated_expressions = {e2, e3, e4}
        
        # Check if this ordering yields the two target expression levels.
        # A tolerance is used for floating-point comparison.
        found_2 = any(math.isclose(expr, 2.0) for expr in calculated_expressions)
        found_5 = any(math.isclose(expr, 5.0) for expr in calculated_expressions)

        if found_2 and found_5:
            # Solution found. Now format and print the results.
            
            # Find which positions (i and j) give 2.0 and 5.0.
            pos_for_2 = 0
            pos_for_5 = 0
            if math.isclose(e2, 2.0): pos_for_2 = 2
            elif math.isclose(e3, 2.0): pos_for_2 = 3
            elif math.isclose(e4, 2.0): pos_for_2 = 4
                
            if math.isclose(e2, 5.0): pos_for_5 = 2
            elif math.isclose(e3, 5.0): pos_for_5 = 3
            elif math.isclose(e4, 5.0): pos_for_5 = 4

            # Create string representations for fractional values for clearer output.
            p_strings = []
            for val in p:
                if math.isclose(val, 1/3): p_strings.append("1/3")
                elif math.isclose(val, 1/2): p_strings.append("1/2")
                elif math.isclose(val, 3/2): p_strings.append("3/2")
                else: p_strings.append(str(int(val)))
            
            print(f"The unique ordering of SNP aFC values is: {', '.join(p_strings)}")
            print("\nThis ordering explains the observed expression levels of the two F2 individuals:")
            
            # Print the equation for the individual with total expression = 2.0
            term1_list = p_strings[pos_for_2 - 1 : 5]
            term2_list = p_strings[0 : pos_for_2]
            print(f"\nIndividual 1 (Total Expression = 2): M/M SNP at position {pos_for_2}")
            print(f"Equation: ( {' * '.join(term1_list)} ) + ( {' * '.join(term2_list)} ) = 2")
            
            # Print the equation for the individual with total expression = 5.0
            term1_list_5 = p_strings[pos_for_5 - 1 : 5]
            term2_list_5 = p_strings[0 : pos_for_5]
            print(f"\nIndividual 2 (Total Expression = 5): M/M SNP at position {pos_for_5}")
            print(f"Equation: ( {' * '.join(term1_list_5)} ) + ( {' * '.join(term2_list_5)} ) = 5")

            # Determine the final answer string based on ranks.
            rank_sequence = [ranks[val] for val in p]
            if rank_sequence[0] < rank_sequence[4]:
                result_str = "".join(map(str, rank_sequence))
            else:
                result_str = "".join(map(str, reversed(rank_sequence)))
            
            print(f"\nThe corresponding sequence of aFC ranks, starting with the lower end rank, is {result_str}.")
            print(f"<<<{result_str}>>>")
            return

# Run the solver
find_snp_order()
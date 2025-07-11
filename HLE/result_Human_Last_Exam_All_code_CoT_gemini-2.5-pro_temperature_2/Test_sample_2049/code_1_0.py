import math
from collections import Counter

def calculate_isomorphism_classes():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This problem is equivalent to counting the number of quintic (degree 5)
    étale algebras over Q that are unramified outside the prime 2.
    An étale algebra of degree 5 is a product of number fields where the sum
    of degrees is 5.

    The first step is to find the number of number fields of degrees 1 to 5
    that are unramified outside of prime 2. These values are obtained from
    number theory databases (like LMFDB) and are as follows:
    - N_1 (degree 1): 1 field (the rational numbers, Q)
    - N_2 (degree 2): 3 fields (Q(i), Q(sqrt(2)), Q(sqrt(-2)))
    - N_3 (degree 3): 0 fields
    - N_4 (degree 4): 11 fields
    - N_5 (degree 5): 6 fields
    """
    
    # N[d] stores the number of fields of degree d unramified outside {2, infinity}.
    # N[0] is a placeholder.
    N = [0, 1, 3, 0, 11, 6]

    # Partitions of the integer 5 represent the structures of the étale algebras.
    partitions = [
        [5],
        [4, 1],
        [3, 2],
        [3, 1, 1],
        [2, 2, 1],
        [2, 1, 1, 1],
        [1, 1, 1, 1, 1]
    ]

    total_count = 0
    calculation_parts = []

    print("The total number is the sum of counts for each partition of 5:\n")

    for p in partitions:
        # Count the multiplicity of each degree in the partition
        counts = Counter(p)
        
        # The number of ways to form an algebra for a given partition is the product
        # of combinations with repetition for each degree.
        # For a degree d with multiplicity m, the number of choices is C(N[d] + m - 1, m).
        term_count = 1
        for degree, multiplicity in counts.items():
            num_fields = N[degree]
            term_count *= math.comb(num_fields + multiplicity - 1, multiplicity)
            
        calculation_parts.append(str(term_count))
        total_count += term_count

        # Explanation for each partition
        partition_str = " + ".join(map(str, p))
        explanation = f"Partition [{partition_str}]: "
        
        formula_parts = []
        for degree, multiplicity in sorted(counts.items()):
            formula_parts.append(f"C({N[degree]} + {multiplicity} - 1, {multiplicity})")
        
        explanation += " * ".join(formula_parts) + f" = {term_count}"
        print(explanation)

    # Final summary of the calculation
    print("\n" + " + ".join(calculation_parts) + f" = {total_count}")

calculate_isomorphism_classes()
<<<27>>>
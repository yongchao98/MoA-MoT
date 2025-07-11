import math
from itertools import combinations

def solve():
    """
    This function finds the number of set tuples (S1, S2, S3, S4) satisfying
    S1 < S2 < S3 < S4 < {1,2,3,4,5} and i in Si for i=1,2,3, where '<' denotes a strict subset.
    
    The problem is solved by distributing the elements {1,2,3,4,5} into 5 regions
    defined by the set chain (R1=S1, R2=S2\S1, etc.) and using the
    Principle of Inclusion-Exclusion to handle the strict inclusion constraints.
    """

    def calculate_assignments(forbidden_regions):
        """
        Calculates the number of ways to assign elements {1,..,5} to regions {R1,..,R5}
        given a set of forbidden regions. Regions are 0-indexed (R1=0, R2=1, ...).
        """
        # Allowed regions for each element {1, 2, 3, 4, 5} based on the problem constraints.
        allowed_sets = [
            {0},                # Element 1 must be in R1
            {0, 1},             # Element 2 can be in R1 or R2
            {0, 1, 2},          # Element 3 can be in R1, R2, or R3
            {0, 1, 2, 3, 4},    # Element 4 can be in any region
            {0, 1, 2, 3, 4}     # Element 5 can be in any region
        ]
        
        num_assignments = 1
        for i in range(5):
            # The number of choices for element i+1 is the size of its allowed set
            # minus any regions that are forbidden in this calculation.
            choices = len(allowed_sets[i] - forbidden_regions)
            num_assignments *= choices
        return num_assignments

    # The strict inclusion requires regions R2, R3, R4, R5 to be non-empty.
    # In our 0-indexed system, these are regions 1, 2, 3, 4.
    non_empty_regions_indices = [1, 2, 3, 4]

    # T: Total assignments satisfying element constraints, ignoring non-empty conditions.
    total_assignments = calculate_assignments(set())

    # S1: Sum of cases where exactly one specified region is empty.
    s1_sum = 0
    s1_terms = []
    for i in non_empty_regions_indices:
        term = calculate_assignments({i})
        s1_terms.append(term)
        s1_sum += term

    # S2: Sum of cases where exactly two specified regions are empty.
    s2_sum = 0
    s2_terms = []
    for i, j in combinations(non_empty_regions_indices, 2):
        term = calculate_assignments({i, j})
        s2_terms.append(term)
        s2_sum += term

    # S3: Sum of cases where exactly three specified regions are empty.
    s3_sum = 0
    s3_terms = []
    for i, j, k in combinations(non_empty_regions_indices, 3):
        term = calculate_assignments({i, j, k})
        s3_terms.append(term)
        s3_sum += term

    # S4: Case where all four specified regions are empty.
    s4_term = calculate_assignments(set(non_empty_regions_indices))

    # Final result from the Principle of Inclusion-Exclusion.
    result = total_assignments - s1_sum + s2_sum - s3_sum + s4_term

    # Format the terms for printing the full equation.
    s1_str = " + ".join(map(str, s1_terms))
    s2_str = " + ".join(map(str, s2_terms))
    s3_str = " + ".join(map(str, s3_terms))
    
    print("The number of valid set arrangements is calculated using the Principle of Inclusion-Exclusion:")
    print(f"Total = T - S1 + S2 - S3 + S4")
    print(f"Total = {total_assignments} - ({s1_str}) + ({s2_str}) - ({s3_str}) + {s4_term}")
    print(f"Total = {total_assignments} - {s1_sum} + {s2_sum} - {s3_sum} + {s4_term}")
    print(f"Total = {result}")

solve()
import math

def solve_number_theory_problem():
    """
    Solves the described number theory problem by counting the number of specified rings
    that have unique length factorizations (are Half-Factorial Domains, HFDs).
    """

    # The problem asks for the size of a set of rings that are HFDs.
    # The set is the union of:
    # A) Rings of integers O_K for K = Q(sqrt(-d)), for d > 0 square-free.
    # B) Rings Z[sqrt(-d)] that are not integrally closed (i.e., d = 3 mod 4).
    # These two sets of rings are disjoint.

    # --- Part 1: Count HFDs in Set A (maximal orders) ---
    # A ring of integers O_K is an HFD if and only if its class number h(-d) is 1 or 2.

    # List of d > 0 where h(-d) = 1 (9 values)
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    # List of d > 0 where h(-d) = 2 (18 values)
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    count_maximal_hfd = len(d_h1) + len(d_h2)

    # --- Part 2: Count HFDs in Set B (non-maximal orders Z[sqrt(-d)]) ---
    # We consider R = Z[sqrt(-d)] where d is square-free and d = 3 (mod 4).
    # Such a ring R is an HFD if and only if its class number, h(R), is <= 2.
    # The class number h(R) can be calculated from the field class number h(-d).

    def get_order_class_number(d, h_K):
        """
        Calculates the class number of the order R = Z[sqrt(-d)] for d = 3 (mod 4).
        h_K is the class number of the maximal order O_K.
        """
        # Special case for d=3. The units of the maximal order are different.
        if d == 3:
            # For d=3, R=Z[sqrt(-3)], h(R)=1.
            return 1
        
        # For d > 3 and d = 3 (mod 4), the formula for h(R) simplifies.
        # The conductor f is 2. The unit index mu_K is 1.
        # h(R) = h_K * (1 - ((-d)/2) / 2) * 2
        
        # We need the Kronecker symbol (-d/2).
        # If d = 3 (mod 8), -d = 5 (mod 8), symbol is -1.
        # If d = 7 (mod 8), -d = 1 (mod 8), symbol is 1.
        kronecker_symbol = 1 if d % 8 == 7 else -1
        
        h_R = h_K * (1 - kronecker_symbol / 2.0) * 2.0
        return int(h_R)

    hfd_non_maximal_d_list = []
    
    # We only need to check d values where h(-d) is small, since h(R) is a multiple of h(-d).
    # If h(-d) > 2, then h(R) > 2. So we check the known lists d_h1 and d_h2.
    all_d_to_check = sorted(d_h1 + d_h2)
    
    for d in all_d_to_check:
        if d % 4 == 3: # We are in the case for non-maximal orders
            h_K = 1 if d in d_h1 else 2
            h_R = get_order_class_number(d, h_K)
            if h_R <= 2:
                hfd_non_maximal_d_list.append(d)

    count_non_maximal_hfd = len(hfd_non_maximal_d_list)

    # --- Final Calculation ---
    total_count = count_maximal_hfd + count_non_maximal_hfd

    print("The total number of rings with unique length factorizations is the sum of:")
    print("1. The number of maximal orders (rings of integers) that are HFDs.")
    print("2. The number of specified non-maximal orders that are HFDs.")
    print("\nCalculation:")
    print(f"Number of HFD maximal orders = {len(d_h1)} (for h=1) + {len(d_h2)} (for h=2) = {count_maximal_hfd}")
    print(f"Number of HFD non-maximal orders (Z[sqrt(-d)] for d=3,7,15) = {count_non_maximal_hfd}")
    print("\nFinal Equation:")
    print(f"{count_maximal_hfd} + {count_non_maximal_hfd} = {total_count}")


solve_number_theory_problem()
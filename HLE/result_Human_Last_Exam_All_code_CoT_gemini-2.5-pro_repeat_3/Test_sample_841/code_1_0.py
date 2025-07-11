import sys

def solve_number_theory_problem():
    """
    This function solves the problem by using known results from algebraic number theory
    regarding the class numbers of imaginary quadratic fields.

    The problem asks for the size of a set of rings that are Half-Factorial Domains (HFDs),
    meaning they have unique lengths of factorizations into irreducibles.

    The set of rings is the union of:
    1. Maximal orders O(Q(sqrt(-d))) for square-free d > 0. These are HFDs if the class number h(-d) <= 2.
    2. Non-maximal orders Z[sqrt(-d)] for square-free d > 0 with d = 3 (mod 4). These are HFDs if h(-d) = 2.

    The two sets of rings are disjoint. We count the HFDs in each set and sum the results.
    """

    # List of square-free d > 0 where the class number h(Q(sqrt(-d))) = 1.
    # These correspond to Unique Factorization Domains (UFDs), which are a subset of HFDs.
    # Source: Baker-Heegner-Stark theorem.
    d_class_number_1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]

    # List of square-free d > 0 where the class number h(Q(sqrt(-d))) = 2.
    # Source: Complete list determined by Baker and Stark.
    d_class_number_2 = [
        5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427
    ]

    # --- Calculation for Maximal Orders ---
    # The maximal orders O(Q(sqrt(-d))) are HFDs if h(-d) is 1 or 2.
    # The number of such rings is the sum of the lengths of the two lists.
    num_maximal_hfd = len(d_class_number_1) + len(d_class_number_2)
    
    # --- Calculation for Non-Maximal Orders ---
    # The non-maximal orders Z[sqrt(-d)] in our set are those with d = 3 (mod 4).
    # They are HFDs if h(-d) = 2.
    # We filter the d_class_number_2 list for this condition.
    non_maximal_hfd_d_values = [d for d in d_class_number_2 if d % 4 == 3]
    num_non_maximal_hfd = len(non_maximal_hfd_d_values)

    # --- Total Count ---
    # The total size is the sum of the counts from the two disjoint sets of rings.
    total_size = num_maximal_hfd + num_non_maximal_hfd
    
    # --- Outputting the Result ---
    # The prompt requires printing the final equation with each number.
    print("This problem requires counting two types of rings from the specified set that have unique factorization lengths (are Half-Factorial Domains).")
    print("\n1. Maximal integer rings O(Q(sqrt(-d))):")
    print(f"   - These are HFDs if the class number is 1 or 2.")
    print(f"   - Number of rings with class number 1: {len(d_class_number_1)}")
    print(f"   - Number of rings with class number 2: {len(d_class_number_2)}")
    print(f"   - Total maximal HFD rings: {len(d_class_number_1)} + {len(d_class_number_2)} = {num_maximal_hfd}")

    print("\n2. Non-maximal integer rings Z[sqrt(-d)] (for d = 3 mod 4):")
    print(f"   - These are HFDs if the class number is exactly 2.")
    print(f"   - Number of rings with class number 2 and d = 3 mod 4: {num_non_maximal_hfd}")
    
    print("\nTotal size is the sum of these disjoint sets.")
    print(f"Final Equation: {num_maximal_hfd} + {num_non_maximal_hfd} = {total_size}")
    
    # Returning the final answer in the specified format
    # The script should not prompt for input, so we use a placeholder or known value.
    # In this context, the value is the calculated total_size.
    # We use file=sys.stdout to ensure it prints to the standard output.
    print(f"\n<<<{total_size}>>>", file=sys.stdout)

solve_number_theory_problem()
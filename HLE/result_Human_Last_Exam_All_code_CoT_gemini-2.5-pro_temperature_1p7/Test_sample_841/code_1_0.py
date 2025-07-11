import math

def main():
    """
    Calculates the size of the specified set of rings with unique factorization lengths.
    """
    
    # Step 1: Define the lists of d for which the maximal order O(Q(sqrt(-d)))
    # has class number (h) equal to 1 or 2. These rings are Half-Factorial Domains (HFDs).
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    # These are the maximal orders in our set.
    count_maximal_hfd = len(d_h1) + len(d_h2)
    
    print("Part 1: Maximal Orders (Rings of Integers)")
    print(f"Found {len(d_h1)} maximal orders with class number 1 for d in {d_h1}.")
    print(f"Found {len(d_h2)} maximal orders with class number 2 for d in {d_h2}.")
    print(f"Total maximal orders with unique factorization lengths: {count_maximal_hfd}\n")
    
    # Step 2: Identify the non-maximal orders Z[sqrt(-d)] that are HFDs.
    # These only occur when d = 3 (mod 4).
    # We analyze them based on d mod 8.
    
    non_maximal_hfd_d = []
    
    # Case d=3: The class number of Z[sqrt(-3)] is 1.
    non_maximal_hfd_d.append(3)
    
    # Case d > 3, d = 3 (mod 8): The class number is 3*h(-d), which is >= 3. No HFDs here.
    
    # Case d = 7 (mod 8): The class number of Z[sqrt(-d)] equals h(-d).
    # We need to find d from our lists where d = 7 (mod 8).
    
    # Combine d values from h=1 and h=2 lists.
    all_d_for_maximal_hfd = d_h1 + d_h2
    
    for d in all_d_for_maximal_hfd:
        # Check for d=7(mod 8), excluding d=3 which was a special case.
        if d != 3 and d % 8 == 7:
            # The class number h(Z[sqrt(-d)]) is h(-d), which is 1 or 2 by definition.
            # So these are HFDs.
            non_maximal_hfd_d.append(d)

    count_non_maximal_hfd = len(non_maximal_hfd_d)

    print("Part 2: Non-Maximal Orders (Z[sqrt(-d)] not integrally closed)")
    print(f"Found {count_non_maximal_hfd} non-maximal orders with unique factorization lengths for d in {sorted(non_maximal_hfd_d)}.")
    
    # Step 3: Sum the counts. The sets of maximal and non-maximal orders are disjoint.
    total_count = count_maximal_hfd + count_non_maximal_hfd
    
    print("\n--- Total Calculation ---")
    print(f"The number of qualifying maximal orders is {count_maximal_hfd}.")
    print(f"The number of qualifying non-maximal orders is {count_non_maximal_hfd}.")
    print(f"The total size of the subset is the sum of these two disjoint sets.")
    print(f"Final equation: {count_maximal_hfd} + {count_non_maximal_hfd} = {total_count}")

if __name__ == "__main__":
    main()
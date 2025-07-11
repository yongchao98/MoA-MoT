def solve_ring_problem():
    """
    This function calculates the size of the set of quadratic integer rings
    (as described in the problem) that are half-factorial domains (HFDs).
    """

    # d values for which the class number of Q(sqrt(-d)) is 1.
    h1_d_values = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    
    # d values for which the class number of Q(sqrt(-d)) is 2.
    h2_d_values = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
    
    # Combined list for d where class number is <= 2.
    h_le_2_d_values = h1_d_values + h2_d_values

    # Count 1: Maximal orders Z[sqrt(-d)] for d = 1, 2 (mod 4).
    # These are HFDs if the class number is <= 2.
    count1 = 0
    for d in h_le_2_d_values:
        if d % 4 == 1 or d % 4 == 2:
            count1 += 1

    # Count 2: Maximal orders Z[(1+sqrt(-d))/2] for d = 3 (mod 4).
    # These are HFDs if the class number is <= 2.
    count2 = 0
    for d in h_le_2_d_values:
        if d % 4 == 3:
            count2 += 1

    # Count 3: Non-maximal orders Z[sqrt(-d)] for d = 3 (mod 4).
    # These are HFDs if (d=3) or (d=7 (mod 8) and class number <= 2).
    # We use a set to avoid double-counting if a value satisfies multiple conditions.
    d_for_s3 = set()
    
    # Condition: d=3. h(-3)=1, so 3 is in our list.
    if 3 in h_le_2_d_values:
      d_for_s3.add(3)
      
    # Condition: d = 7 (mod 8) and class number <= 2.
    for d in h_le_2_d_values:
        if d % 8 == 7:
            d_for_s3.add(d)
            
    count3 = len(d_for_s3)

    # The total size is the sum of the counts from the three disjoint sets.
    total_count = count1 + count2 + count3
    
    # Print the breakdown of the final sum as requested.
    print(f"{count1} + {count2} + {count3} = {total_count}")

solve_ring_problem()
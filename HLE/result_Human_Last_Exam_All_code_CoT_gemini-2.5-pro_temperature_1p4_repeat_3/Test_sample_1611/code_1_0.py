import collections

def solve():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid.
    """
    n = 4
    m = 4

    # Step 1: Find all valid row patterns for a row of length m.
    # A pattern is valid if no two adjacent bits are 1.
    # We can check this with the bitwise operation (p & (p << 1)) == 0.
    valid_patterns = []
    for p in range(1 << m):
        if (p & (p << 1)) == 0:
            valid_patterns.append(p)

    num_patterns = len(valid_patterns)
    pattern_map = {pattern: i for i, pattern in enumerate(valid_patterns)}

    print(f"For a {n}x{m} grid.")
    print(f"There are {num_patterns} valid patterns for a single row:")
    # Printing in binary format for clarity
    print([f"{p:0{m}b}" for p in valid_patterns])
    print("-" * 30)

    # Step 2: Create a compatibility map.
    # Two patterns are compatible if they can be placed in adjacent rows.
    # This means their bitwise AND is 0.
    compatibility = collections.defaultdict(list)
    for p1 in valid_patterns:
        idx1 = pattern_map[p1]
        for p2 in valid_patterns:
            idx2 = pattern_map[p2]
            if (p1 & p2) == 0:
                compatibility[idx1].append(idx2)

    # Step 3: Dynamic Programming calculation.
    # dp[i] will store the number of ways to fill a grid of i rows ending with pattern valid_patterns[i]
    
    # Base case: For the first row, there is one way for each valid pattern.
    dp = [1] * num_patterns
    print(f"Row 1: Counts for each pattern = {dp}")
    print(f"Total for 1 row = {sum(dp)}")
    print("-" * 30)

    # Iterate from the second row to the n-th row.
    for i in range(2, n + 1):
        new_dp = [0] * num_patterns
        for j in range(num_patterns):
            # Sum the counts of compatible patterns from the previous row
            for compatible_idx in compatibility[j]:
                new_dp[j] += dp[compatible_idx]
        dp = new_dp
        print(f"Row {i}: Counts for each pattern = {dp}")
        print(f"Total for {i} rows = {sum(dp)}")
        print("-" * 30)

    # Step 4: The final answer is the sum of counts for all patterns in the last row.
    total_assignments = sum(dp)

    # Print the final calculation as an equation
    final_sum_str = " + ".join(map(str, dp))
    print(f"Final calculation for the {n}x{m} grid:")
    print(f"{final_sum_str} = {total_assignments}")
    print("\nSo, the total number of different 0/1 assignments is:")
    print(total_assignments)


solve()
<<<1289>>>
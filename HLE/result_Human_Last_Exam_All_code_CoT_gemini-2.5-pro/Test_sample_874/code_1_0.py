def solve():
    """
    Solves the problem by identifying the optimal tuple based on the Tribonacci sequence
    and then performing the required calculation.
    """
    # Step 2 & 3: Generate Tribonacci numbers up to the limit 10,000,000
    limit = 10_000_000
    trib = [0, 0, 1]
    while trib[-1] <= limit:
        next_trib = trib[-1] + trib[-2] + trib[-3]
        if next_trib > limit:
            break
        trib.append(next_trib)

    # The optimal tuple uses the three largest available Tribonacci numbers and 0.
    # The list 'trib' contains T_0, T_1, ..., T_31.
    # The specific arrangement (0, T_n, T_{n+1}, T_{n+2}) maximizes the sequence length.
    # Here, T_{n+2} must be the largest, so T_{n+2} = T_31. This means n = 29.
    # Our tuple's values come from {0, T_29, T_30, T_31}.
    n = len(trib) - 3  # This corresponds to index 29 for T_29

    # Step 4: Identify the specific tuple (a, b, c, d).
    # We choose the lexicographically smallest tuple among those that give the maximum length.
    # This corresponds to the canonical form (0, T_n, T_{n+1}, T_{n+2}).
    a = 0
    b = trib[n]      # T_29
    c = trib[n + 1]  # T_30
    d = trib[n + 2]  # T_31

    # Step 5: Compute the expression (a + b - c - d) mod 1000.
    # Using integer arithmetic to avoid potential floating point issues.
    result = (a + b - c - d) % 1000

    # The problem asks to output the numbers in the final equation.
    print(f"{a} + {b} - {c} - {d} mod 1000 = {result}")
    
    # Also printing the final answer in the requested format.
    print(f"\n<<<{result}>>>")

solve()

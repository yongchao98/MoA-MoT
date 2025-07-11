def solve():
    """
    This function calculates the number of partitions of an integer n.
    This corresponds to the number of conjugacy classes of the symmetric group S_n,
    which is the answer to the problem for n.
    """
    n = 9

    # dp[i] will store the number of partitions of integer i.
    # We initialize dp array of size n+1 with all zeros.
    dp = [0] * (n + 1)
    
    # There is one way to partition 0 (the empty sum).
    dp[0] = 1

    # We build the partition table using dynamic programming.
    # We iterate through all numbers from 1 to n to be used as parts.
    for i in range(1, n + 1):
        # For each part 'i', we update the partition counts for numbers j >= i.
        # This is because 'i' can be a part of a partition of any number j >= i.
        # dp[j] is updated by adding dp[j-i], which represents the number of ways
        # to form the rest of the partition.
        for j in range(i, n + 1):
            dp[j] += dp[j - i]
    
    # The final answer is the number of partitions of n.
    result = dp[n]
    
    # As requested, printing the final result.
    # For this problem, the final "equation" is just the resulting number.
    print(f"The number of partitions of {n} is p({n}).")
    print(f"p(1) = {dp[1]}")
    print(f"p(2) = {dp[2]}")
    print(f"p(3) = {dp[3]}")
    print(f"p(4) = {dp[4]}")
    print(f"p(5) = {dp[5]}")
    print(f"p(6) = {dp[6]}")
    print(f"p(7) = {dp[7]}")
    print(f"p(8) = {dp[8]}")
    print(f"p(9) = {result}")

solve()
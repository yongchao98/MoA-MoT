def count_stable_reductions(g):
    """
    Calculates the number of types of stable reductions of a genus g curve
    whose Jacobian has good reduction. This is equivalent to calculating the
    number of partitions of the integer g.
    
    This function uses a dynamic programming approach.
    dp[i] will store the number of partitions of the integer i.
    """
    
    if g < 0:
        return 0
    
    # dp[i] will be storing the number of partitions of i
    dp = [0] * (g + 1)
    dp[0] = 1 # There is one partition of 0 (the empty sum)
    
    # Build up the partition counts
    # For each number i from 1 to g (representing the parts of the partition)
    for i in range(1, g + 1):
        # Update the dp array for numbers j >= i
        for j in range(i, g + 1):
            dp[j] += dp[j - i]
            
    return dp[g]

# The genus of the curve is 4.
genus = 4
num_types = count_stable_reductions(genus)

# The corresponding partitions represent the distribution of genus
# among the components of the reduced curve.
partitions_str = ["4", "3 + 1", "2 + 2", "2 + 1 + 1", "1 + 1 + 1 + 1"]

print(f"The problem is equivalent to counting the number of integer partitions of the genus, g = {genus}.")
print("The types of stable reductions correspond to the following partitions of the genus:")
for p in partitions_str:
    print(f"- A configuration with component genera summing to 4: {p}")

print(f"\nThe total number of types is the number of partitions of {genus}, which is {num_types}.")

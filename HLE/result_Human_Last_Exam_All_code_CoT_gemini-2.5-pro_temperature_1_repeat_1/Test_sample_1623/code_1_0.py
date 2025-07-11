# Grid size
n = 5

# The maximal Thurston-Bennequin number for the knot type is found by 
# minimizing the number of inversions for a permutation with the same cycle structure.
# The given permutation is a 5-cycle. The minimum number of inversions for a 5-cycle in S_5 is 4.
min_inversions = 4

# Calculate the total number of pairs (i, j) with i < j
num_pairs = n * (n - 1) // 2

# Calculate the maximal writhe using the minimal number of inversions
# writhe = num_ascents - num_inversions = (num_pairs - inv) - inv = num_pairs - 2 * inv
max_writhe = num_pairs - 2 * min_inversions

# Calculate the maximal Thurston-Bennequin number using the formula tb = writhe - n
max_tb = max_writhe - n

# Output the final equation with each number
print(f"{max_writhe} - {n} = {max_tb}")
